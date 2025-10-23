'use strict';

import * as THREE from 'three';
import { ViewManager } from './view.js';

const EARTH_RADIUS_KM  = 6371;
const KM_PER_UNIT = EARTH_RADIUS_KM / 50; // Earth radius to units
const fromUnits = u => u * KM_PER_UNIT;
const toDegrees = rad => rad * 180 / Math.PI;

// Standard gravitational parameter for Earth in km^3/s^2
const MU_EARTH = 398600.4418;

const temp1Vec = new THREE.Vector3();
const temp2Vec = new THREE.Vector3();

// Orient the local coordinate frame of an object on the surface of
// a parent body so that local XZ plane is tangent to surface and
// +X is the local meridian pointing North
function orientToSurface(marker) {
  // Point on the surface must be in parent coordinates
  const p = marker.position;

  // Local UP is the outward normal from sphere
  const y_local = p.clone().normalize();

  // Calculate tangential velocity ω × p with ω = (0, 1, 0) for Y-axis
  // rotation and p(x, y, z) --> (z, 0, -x)
  const v_east = new THREE.Vector3(p.z, 0, -p.x);

  // Detect special case of point exactly on the global Y-axis (pole) and
  // fallback on parent X-axis
  const onY = v_east.lengthSq() === 0;
  const z_local = onY ? new THREE.Vector3(1, 0, 0) : v_east.normalize();

  // North is the cross product up X east vectors
  const x_local = new THREE.Vector3().crossVectors(y_local, z_local);

  // Rotation matrix expresses the local coordinate axes in the parent frame
  // and maps a point's local frame to its parent's frame: vworld​ = R ⋅ vlocal​
  const rotationMatrix = new THREE.Matrix4();
  rotationMatrix.makeBasis(x_local, y_local, z_local);

  // Set point's quaternion from this rotation matrix, rendering will
  // rotate the marker so that its local coordinate axes align with the
  // parent-space directions given by x_local, y_local, and z_local
  marker.quaternion.setFromRotationMatrix(rotationMatrix);
  marker.updateWorldMatrix(true, false);
}

// Orient satellite local frame so camera (that by default
// looks toward z-axis) will point "in the rear-view mirror"
function orientSatellite(satellite, earth) {
  const earthPosWorldVec = new THREE.Vector3();
  earth.getWorldPosition(earthPosWorldVec);
  // Rotates satellite so that its +Z axis points towards earth
  satellite.lookAt(earthPosWorldVec);
  // Now rotate to point to horizon 'behind'
  const R = earth.geometry.parameters.radius;
  const horizonAngle = Math.asin(R / satellite.position.length());
  satellite.rotateY(horizonAngle);
}

// Calculates the Azimuth and Elevation of the views's current target point
// and Latitude and Longitude of observer position.
// The Az+ El are in the local reference frame of the marker object, instead
// Lat + Lon are in geocentric coordinate system.
function getGeoData(marker, targetWorldPos) {

  const targetLocalPos = marker.worldToLocal(temp1Vec.copy(targetWorldPos));

  // If the target is exactly at the observer's position, the direction is undefined.
  if (targetLocalPos.lengthSq() === 0) {
    return [];
  }

  // Elevation is the angle between the target vector and its projection on the XZ plane.
  // We use asin(y / length).
  const elevationRad = Math.asin(targetLocalPos.y / targetLocalPos.length());

  // Azimuth is the angle in the XZ plane, measured _clockwise_ from North (+X).
  // The standard function atan2(y, x) measures a counter-clockwise angle in the
  // XY plane when viewed from the +Z axis.
  // By using atan2(z, x), we are calculating the angle in the XZ plane. A standard
  // right-hand rotation transforms the original +Z viewing axis to our -Y axis.
  // This means the rotation is counter-clockwise when viewed from -Y.
  // Therefore, when we view it from our standard +Y (top-down) direction, the rotation
  // is perceived as CLOCKWISE, which is exactly what a navigational azimuth requires.
  const x_to_z = Math.atan2(targetLocalPos.z, targetLocalPos.x);

  // Normalize if the result is negative (target is in the "western" half of the circle)
  const azimuthRad = x_to_z < 0 ? x_to_z + 2 * Math.PI : x_to_z;

  // The marker's position is already in the local coordinate system of the
  // parent object (the planet), which is our geocentric frame.
  const observerLocalPos = marker.position;

  // Latitude: The angle above or below the body's equatorial (XZ) plane.
  // We use asin(y / length) to get the geocentric latitude.
  const latitudeRad = Math.asin(observerLocalPos.y / observerLocalPos.length());

  // Longitude: The angle around the polar (Y) axis in the equatorial plane.
  // We assume the prime meridian (0° longitude) is along the parent's +X axis.
  // Using atan2(z, x) gives a counter-clockwise angle from +X, matching the
  // convention where East longitudes are positive and West are negative.
  const longitudeRad = Math.atan2(observerLocalPos.z, observerLocalPos.x);

  return {
    azimuth:   toDegrees(azimuthRad),
    elevation: toDegrees(elevationRad),
    latitude:  toDegrees(latitudeRad),
    longitude: toDegrees(longitudeRad),
  };
}

let nadirMarker = null;
let prevElapsed = 0;
const prevNadirLocalPos = new THREE.Vector3();
let prevGroundSpeed = 0;
let prevInertialSpeed = 0;
let prevHeight = 0;

// Satellite position must be in Earth local frame
function placeAtNadir(nadirMarker, radius, curSatPos) {
  const nadirPos = temp1Vec.copy(curSatPos).normalize().multiplyScalar(radius);
  nadirMarker.position.copy(nadirPos);
  return nadirMarker.position;
}

// Calculates satellite ground speed and other surface data
function getSatData(satellite, earth, elapsed) {

  const earthRadius = earth.geometry.parameters.radius;
  const curSatPos = satellite.position;

  // Marker on the ground used for ground speed
  if (!nadirMarker) {
      nadirMarker = new THREE.Object3D();
      earth.add(nadirMarker);
      const curPos = placeAtNadir(nadirMarker, earthRadius, curSatPos);
      prevNadirLocalPos.copy(curPos);
      prevElapsed = elapsed;
      return null;
  }

  // We can move in both time directions
  const deltaTime = Math.abs(elapsed - prevElapsed);
  if (deltaTime === 0) { // paused
    return {
      height: prevHeight,
      inertialSpeed: prevInertialSpeed,
      groundSpeed: prevGroundSpeed,
    }
  }

  // Compute orbital speed by Vis-viva equation
  const r_km = fromUnits(curSatPos.length()); // in km
  const a_km = EARTH_RADIUS_KM + 417; // FIXME do not hardcode
  const speed_km_per_s_squared = MU_EARTH * ((2 / r_km) - (1 / a_km));
  const speed_km_per_s = Math.sqrt(speed_km_per_s_squared);
  const inertialSpeed = Math.floor(speed_km_per_s * 3600);

  // Place nadir marker to new current position
  const curNadirPos = placeAtNadir(nadirMarker, earthRadius, curSatPos);

  // Get the angle between the previous and current nadir position vectors.
  const angle = prevNadirLocalPos.angleTo(curNadirPos);

  // The distance is the arc length: s = r * θ
  const groundDistance = fromUnits(earthRadius) * angle;
  const groundSpeed = Math.floor(3600 * groundDistance / deltaTime); // km/h

  // Compute height above the ground
  const heightVec = temp1Vec.copy(curSatPos).sub(curNadirPos);
  const height = Math.floor(fromUnits(heightVec.length()));

  prevElapsed = elapsed;
  prevNadirLocalPos.copy(curNadirPos);
  prevGroundSpeed = groundSpeed;
  prevInertialSpeed = inertialSpeed;
  prevHeight = height;

  return {
    height: height,
    inertialSpeed: inertialSpeed,
    groundSpeed: groundSpeed,
   };
}

// Tracks state for an observer on a planet or on a satellite
export class Observer {
    constructor(views, onSatellite) {
        this.views = views;
        this.object = null; // usually Earth
        this.marker = null; // in local coordinates
        this.viewIndex = null;
        this.onSatellite = onSatellite;
    }

    getView() {
      return this.viewIndex !== null ? this.views.get(this.viewIndex) : null;
    }

    getGeoData() {
      const view = this.getView();
      return getGeoData(this.marker, view.controls.target);
    }

    getSatData(elapsedSeconds) {
      if (this.object === null)
        return null;
      return getSatData(this.marker, this.object, elapsedSeconds);
    }

    // Set the oberver's local view
    createObserverView(object, marker, eyeHeight, lookAhead) {

      // Currently we handle only one observer
      console.assert(this.object === null, "Setting already existing observer");

      // Marker position should be already in local frame
      if (this.onSatellite){
        orientSatellite(marker, object);
      } else {
        orientToSurface(marker);
      }
      this.object = object;
      this.marker = marker;

      // Create a new observer view placed on the marker
      const fov = this.onSatellite ? 90 : 110;
      const viewIndex = this.views.createObserverView(marker, eyeHeight, lookAhead, fov);
      this.viewIndex = viewIndex;

      // Lock the view camera on the marker
      const view = this.views.get(viewIndex);
      view.lockTo(marker, true, this.onSatellite);
      return viewIndex;
    }

    // Move the marker to a new position on planet surface
    placeAt(latDeg, lonDeg) {

      const view = this.views.get(this.viewIndex);
      if (this.views.getActive() !== view)
        return;

      const latRad = toRadians(latDeg);
      const lonRad = toRadians(lonDeg);
      const radius = this.marker.position.length();

      // Convert spherical coordinates to Cartesian (x, y, z) in the planet's local
      // frame. Y is the polar axis. +X is the prime meridian (0° longitude).
      // The -z for sin(lat) is to match the longitude calculation in view.getGeoData
      const newSurfacePoint = new THREE.Vector3();
      newSurfacePoint.y = radius * Math.sin(latRad);
      const xzRadius = radius * Math.cos(latRad); // Radius of the circle at this latitude
      newSurfacePoint.x = xzRadius * Math.cos(lonRad);
      newSurfacePoint.z = -xzRadius * Math.sin(lonRad);

      // Move the marker to new surface point and realign its plane
      // New marker's position must be in local coordinates
      this.marker.position.copy(newSurfacePoint);
      orientToSurface(this.marker);

      // Camera is already locked to the marker, locking update will
      // set correct camera and target position at next frame.
    }

    // Drop observer view and release relative resources
    disposeObserverView() {
      const view = this.getView();
      view.unlock();
      this.views.dispose(this.viewIndex);
      this.object = null;
      this.viewIndex = null;
      this.marker = null;
    }
};
