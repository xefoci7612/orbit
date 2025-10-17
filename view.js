'use strict';

import * as THREE from 'three';
import { OrbitControls } from 'orbitControls';

// ============================================================================
// CAMERA VIEWS
// ============================================================================

const tempVec = new THREE.Vector3();
const tempQuat = new THREE.Quaternion();

// Helper to convert radians to degrees
const toDegrees = rad => rad * 180 / Math.PI;

class CameraLock {
  constructor(marker, view, isFixedCamera) {
    this.target = marker; // The object on the planet's surface the camera is attached to
    this.view = view;
    this.object = marker.parent;
    this.isFixedCamera = isFixedCamera;
    this.prevPosition = new THREE.Vector3();
    this.prevOrientation = new THREE.Quaternion();
    this.cameraLocalPos = new THREE.Vector3(); // stable with planet rotation

    this.target.getWorldPosition(this.prevPosition);
    this.target.getWorldQuaternion(this.prevOrientation);
    this.cameraLocalPos.copy(this.view.camera.position);
    this.target.worldToLocal(this.cameraLocalPos);
  }

  update() {
    // Camera position is always in world coords
    const currentPos = this.target.getWorldPosition(tempVec);
    const currentQuat = this.target.getWorldQuaternion(tempQuat);

    // Detect target movement since last frame
    const translationDelta = currentPos.clone().sub(this.prevPosition);
    const rotationDelta = currentQuat.clone().multiply(this.prevOrientation.clone().invert());

    // Cache for next frame
    this.prevPosition.copy(currentPos);
    this.prevOrientation.copy(currentQuat);

    // Applies a rigid-body transform delta (translation + rotation) to the active camera
    // while preserving the user's relative viewpoint (offset from target).
    // This is used to keep the camera locked to a moving/rotating object.
    const camera = this.view.camera;
    const controls = this.view.controls;

    // Preserve current offset from target (in world space)
    const targetToCamera = camera.position.clone().sub(controls.target);

    // Rotate the offset vector by the object's rotation delta
    targetToCamera.applyQuaternion(rotationDelta);

    // Move the target by the object's translation delta
    controls.target.add(translationDelta);

    // Reconstruct camera position
    camera.position.copy(controls.target).add(targetToCamera);

    controls.update();

    if (this.isFixedCamera && controls.enabled)
      this.targetCompensation();
  }

  targetCompensation() {

    // Camera and target are always in World coordinates
    const camera = this.view.camera;
    const controls = this.view.controls;
    const stableCameraWorldPos = tempVec.copy(this.cameraLocalPos).clone();
    this.target.localToWorld(stableCameraWorldPos);

    // P1: The CORRECT, stable camera position on the planet's surface.
    const P1 = stableCameraWorldPos;

    // P2: The ACTUAL camera position after OrbitControls has moved it.
    const P2 = camera.position;

    // T1: The current target (pivot point) of the OrbitControls.
    const T1 = controls.target;

    if (P1.equals(P2))
      return;

    // ========================================================================
    // PART 1: Compute the user's intended rotation (`rotP1_P2`)
    // We calculate the rotation that OrbitControls applied to the camera
    // by moving it from P1 to P2 around the central target T1.
    // ========================================================================

    // Get the vector from the target to the camera's correct starting position (P1).
    const vec_T1_P1 = tempVec.subVectors(P1, T1).clone();

    // Get the vector from the target to the camera's actual final position (P2).
    const vec_T1_P2 = tempVec.subVectors(P2, T1).clone();

    // We only care about the change in direction, so we normalize the vectors.
    vec_T1_P1.normalize();
    vec_T1_P2.normalize();

    // Calculate the quaternion representing the rotation from the "before" vector
    // to the "after" vector. This is the pure rotation intended by the user.
    const userRotationDelta = tempQuat.setFromUnitVectors(vec_T1_P1, vec_T1_P2).clone();

    // ========================================================================
    // PART 2 & 3: Apply the correction.
    // We force the camera back to P1 and move the target to a new position (T2)
    // to replicate the user's intended rotation from the correct viewpoint.
    // ========================================================================

    // First, calculate the original line-of-sight vector from the correct
    // camera position (P1) to the target (T1).
    const sightVector = tempVec.subVectors(T1, P1).clone();

    // Apply the captured user rotation directly to this line of sight.
    // This gives us the new direction the camera should be pointing.
    sightVector.applyQuaternion(userRotationDelta);

    // The new target (T2) is found by starting at the fixed camera position (P1)
    // and looking out along the new, rotated sight vector.
    const newTargetPosition = tempVec.addVectors(P1, sightVector).clone();

    // --- EXECUTE THE CORRECTION ---
    // Force the camera back to its stable, grounded position.
    camera.position.copy(P1);

    // Up direction drifts with orbit, realign it
    const cameraUpWorld = tempVec.set(0, 1, 0); // Y-axis
    cameraUpWorld.transformDirection(this.target.matrixWorld);
    camera.up.copy(cameraUpWorld);

    // Update the controls' target to the new, compensated position.
    controls.target.copy(newTargetPosition);

    //controls.update(); //no update here!
  }
};

class View {
  constructor(camera, controls, earthPosWorldVec) {
    this.camera = camera;
    this.controls = controls;
    // Initial position and target are geocentric
    this.initialPosition = camera.position.clone().sub(earthPosWorldVec);
    this.initialTarget = controls.target.clone().sub(earthPosWorldVec);
    this.cameraLock = null;
  }

  lockTo(marker, isFixedCamera) {
    // Set a lock on the view so that camera will follow the marker
    this.cameraLock = new CameraLock(marker, this, isFixedCamera);
  }

  unlock() {
    this.cameraLock = null;
  }

  // Calculates the Azimuth and Elevation of the views's current target point
  // and Latitude and Longitude of observer position.
  // The Az+ El are in the local reference frame of the marker object, instead
  // Lat + Lon are in geocentric coordinate system.
  getGeoData(marker) {

    const targetWorldPos = this.controls.target;
    const targetLocalPos = marker.worldToLocal(tempVec.copy(targetWorldPos));

    // If the target is exactly at the observer's position, the direction is undefined.
    if (targetLocalPos.lengthSq() === 0) {
      return [];
    }

    // Elevation is the angle between the target vector and its projection on the XZ plane.
    // We use asin(y / length).
    const elevationRad = Math.asin(targetLocalPos.y / targetLocalPos.length());

    // Azimuth is the angle in the XZ plane.
    // We use atan2(x, z) to get the angle from the +Z axis, growing in
    // anticlockwise direction.
    const azimuthRadMath = Math.atan2(targetLocalPos.x, targetLocalPos.z);

    // Geographical Azimuth instead starts from the -Z axis (North),
    // growing in a clockwise direction.
    const azimuthRad = Math.PI - azimuthRadMath;

    // The marker's position is already in the local coordinate system of the
    // parent object (the planet), which is our geocentric frame.
    const observerLocalPos = marker.position;

    // Latitude: The angle above or below the body's equatorial (XZ) plane.
    // We use asin(y / radius).
    const latitudeRad = Math.asin(observerLocalPos.y / observerLocalPos.length());

    // Longitude: The angle around the polar (Y) axis in the equatorial plane.
    // We assume the prime meridian (0° longitude) is along the parent's +X axis.
    // Using atan2(-z, x) gives a clockwise angle from +X, matching the convention
    // where East longitudes are positive.
    const longitudeRad = Math.atan2(-observerLocalPos.z, observerLocalPos.x);

    return {
      azEl: [toDegrees(azimuthRad), toDegrees(elevationRad)],
      latLon: [toDegrees(latitudeRad), toDegrees(longitudeRad)]
    };
  }
};

class ViewManager {
  constructor(scene, renderer, defaultPosition, earth) {
    this.activeIdx = null;
    this.views = [];
    this.scene = scene;
    this.renderer = renderer;
    this.defaultPosition = defaultPosition; // relative to Earth
    this.earth = earth;
  }

  init(earthPosWorldVec) {
    const target = earthPosWorldVec.clone();
    const position = earthPosWorldVec.clone().add(this.defaultPosition);
    const mainViewIndex = this.create({ position: position,
                                        up: new THREE.Vector3(0, 1, 0) },
                                      { target: target });
    this.setActive(mainViewIndex);
  }

  get(idx) {
    return this.views[idx];
  }

  getActive() {
    return this.activeIdx !== null ? this.views[this.activeIdx] : null;
  }

  getOrbitLockedObjects() {
    return this.views.filter(v => v.cameraLock !== null && !v.cameraLock.isFixedCamera)
                     .map(v => v.cameraLock.object); // ignore observer views
  }

  create(cameraConfig, controlsConfig) {
    const camera = new THREE.PerspectiveCamera(
      cameraConfig.fov || 45,
      innerWidth / innerHeight,
      cameraConfig.near || 0.001,
      cameraConfig.far || 10000
    );

    // orbitControls camera should always be in world coordinates
    this.scene.add(camera);

    // Camera position and target must be in world coordinates
    const controls = new OrbitControls(camera, this.renderer.domElement);
    controls.target.copy(controlsConfig.target);
    camera.position.copy(cameraConfig.position);
    camera.up.copy(cameraConfig.up);
    controls.update(); // Sync controls with initial state
    controls.enabled = false; // New view starts disabled

    // Apply specific properties if they exist in the config
    if (controlsConfig.enablePan !== undefined) controls.enablePan = controlsConfig.enablePan;
    if (controlsConfig.enableZoom !== undefined) controls.enableZoom = controlsConfig.enableZoom;
    if (controlsConfig.minPolarAngle !== undefined) controls.minPolarAngle = controlsConfig.minPolarAngle;
    if (controlsConfig.maxPolarAngle !== undefined) controls.maxPolarAngle = controlsConfig.maxPolarAngle;

    const earthPosWorldVec = new THREE.Vector3();
    this.earth.getWorldPosition(earthPosWorldVec);
    const view = new View(camera, controls, earthPosWorldVec);
    const viewIndex = this.views.push(view) - 1;
    return viewIndex;
  }

  dispose(viewIndex) {
    if (viewIndex === this.activeIdx) {
      this.setDefault(); // dropping view under our feet
    }
    const view = this.get(viewIndex);
    view.camera.parent.remove(view.camera);
    view.controls.dispose();
    console.assert(view.cameraLock === null, "Disposing view with an active lock");
  }

  setActive(viewIndex) {

    // Disable current view
    const curView = this.getActive();
    if (curView) // Can be null at startup
      curView.controls.enabled = false;

    // Switch to active view
    this.activeIdx = viewIndex;

    // Set new view default start position and update controls
    const newView = this.getActive();
    newView.camera.position.copy(newView.initialPosition);
    newView.controls.target.copy(newView.initialTarget);

    // Initial positions are geocentric, convert to world coordinates
    const earthPosWorldVec = new THREE.Vector3();
    this.earth.getWorldPosition(earthPosWorldVec);
    newView.camera.position.add(earthPosWorldVec);
    newView.controls.target.add(earthPosWorldVec);

    newView.controls.enabled = true;
    newView.controls.update();

    // Important: Update camera aspect ratio on switch
    newView.camera.aspect = innerWidth / innerHeight;
    newView.camera.updateProjectionMatrix();
  }

  setDefault() {
    this.setActive(0); // initial view is our default
  }

  clone(view) {
    const camera = view.camera;
    const controls = view.controls;
    const newViewIndex = this.create({ position: camera.position, up: camera.up },
                                     { target: controls.target });
    return newViewIndex;
  }

  update() {
    // Update positions of locked cameras
    this.views.forEach(view => {
      if (view.cameraLock)
        view.cameraLock.update();
    });
    // Update controls of active camera
    const activeView = this.getActive();
    activeView.controls.update();
    return activeView.camera;
  }

  // Lock active camera into a geostationary orbit around the selected object
  lockToOrbitView(marker, objectToLock) {

    const view = this.getActive();
    const camera = view.camera;
    const controls = view.controls;

    // Set camera target on object center for better UX (user can still pan after)
    const temp = new THREE.Vector3();
    const surfacePoint = marker.getWorldPosition(temp).clone();
    const objCenter = objectToLock.getWorldPosition(temp).clone();
    controls.target.copy(objCenter);

    // Reposition camera along center-surface vector, same distance
    const distanceToCenter = camera.position.distanceTo(objCenter);
    const direction = temp.subVectors(surfacePoint, objCenter).normalize();
    camera.position.copy(direction.multiplyScalar(distanceToCenter)).add(objCenter);
    controls.update();
    view.lockTo(marker, false);
  }

  // Setup the camera view for a observer on a planet surface
  createObserverView(marker, eyeHeight, lookAhead) {

    // Camera postion and target must be in world coordinates
    const cameraPosWorld = new THREE.Vector3(0, eyeHeight, 0);
    const targetPosWorld = new THREE.Vector3(0, 5 * eyeHeight, lookAhead);
    marker.localToWorld(cameraPosWorld);
    marker.localToWorld(targetPosWorld);

    // Convert UP direction to world space (rotation only!)
    const cameraUpWorld = new THREE.Vector3(0, 1, 0); // Y-axis
    cameraUpWorld.transformDirection(marker.matrixWorld);

    const cameraConfig = {
        position: cameraPosWorld,
        fov: 60,
        up: cameraUpWorld,
        //rollAngle: Math.PI / 4
    };
    const controlsConfig = {
        target: targetPosWorld,
        enablePan: false,
        enableZoom: false,
        //maxPolarAngle: Math.PI / 2,  // Prevent looking down through the ground
    }

    // Create a new view and set a lock on so that camera will not drift
    const viewIndex = this.create(cameraConfig, controlsConfig);
    return viewIndex;
  }
};

export { ViewManager };
