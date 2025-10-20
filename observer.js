'use strict';

import * as THREE from 'three';
import { ViewManager } from './view.js';

// Orient the local coordinate frame of an object on the surface of
// a parent body so that local XZ plane is tangent to surface with
// X axis in the same direction of planet rotation (toward East)
function orientToSurface(marker) {
  // Point on the surface must be in parent coordinates
  const p = marker.position;

  // Calculate Local Y-axis (outward normal from sphere)
  const y_local = p.clone().normalize();

  // Calculate Local X-axis as tangential velocity (z, 0, -x) for Y-axis rotation
  // Here v⋅p = 0 for both p = Y axis and any radial vector (x, y, z)
  const v = new THREE.Vector3(p.z, 0, -p.x);

  // Detect special case of point exactly on the global Y-axis (pole) and
  // fallback on parent X-axis
  const onY = v.lengthSq() === 0;
  const x_local = onY ? new THREE.Vector3(1, 0, 0) : v.normalize();

  // Local Z-axis is cross product of x_local and y_local for right-handed system
  const z_local = new THREE.Vector3().crossVectors(x_local, y_local);

  // Rotation matrix whose columns are the local basis vectors (x_local, y_local, z_local)
  // This matrix expresses the local coordinate axes in the parent (world) frame
  // and maps a point's local frame to its parent's frame: vworld​ = R ⋅ vlocal​
  const rotationMatrix = new THREE.Matrix4();
  rotationMatrix.makeBasis(x_local, y_local, z_local);

  // Set point's quaternion from this rotation matrix, rendering will
  // rotate the marker so that its local coordinate axes align with the
  // parent-space directions given by x_local, y_local, and z_local
  marker.quaternion.setFromRotationMatrix(rotationMatrix);
  marker.updateWorldMatrix(true, false);
}

// Tracks state for an observer on a planet
export class Observer {
    constructor(views) {
        this.views = views;
        this.object = null;
        this.marker = null; // in local coordinates
        this.viewIndex = null;
        this.tempVec = new THREE.Vector3();
    }

    getView() {
      return this.viewIndex !== null ? this.views.get(this.viewIndex) : null;
    }

    // Set the oberver's local view
    enterObserverView(object, marker, eyeHeight, lookAhead) {

      // Currently we handle only one observer
      console.assert(this.object === null, "Setting already existing observer");

      // Marker position should be already in local frame
      orientToSurface(marker);
      this.object = object;
      this.marker = marker;

      // Create a new observer view placed on the marker
      const viewIndex = this.views.createObserverView(marker, eyeHeight, lookAhead);
      this.viewIndex = viewIndex;

      // Lock the view camera on the marker
      const view = this.views.get(viewIndex);
      view.lockTo(marker, true);
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
    exitObserverView() {
      const view = this.getView();
      view.unlock();
      this.views.dispose(this.viewIndex);
      this.object = null;
      this.viewIndex = null;
      this.marker = null;
    }
};
