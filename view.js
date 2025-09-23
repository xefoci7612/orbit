'use strict';

import * as THREE from 'three';
import { OrbitControls } from 'orbitControls';

// ============================================================================
// CAMERA VIEWS
// ============================================================================

class Lock {
  constructor(marker, view) {
    this.target = marker;
    this.view = view;
    this.object = marker.parent;
    this.prevPosition = new THREE.Vector3();
    this.prevOrientation = new THREE.Quaternion();
    this.tempVec = new THREE.Vector3();
    this.tempQuat = new THREE.Quaternion();

    this.target.getWorldPosition(this.prevPosition);
    this.target.getWorldQuaternion(this.prevOrientation);
  }

  update() {
    // Camera position is always in world coords
    const currentPos = this.target.getWorldPosition(this.tempVec);
    const currentQuat = this.target.getWorldQuaternion(this.tempQuat);

    // Detect target movement since last frame
    const translationDelta = currentPos.clone().sub(this.prevPosition);
    const rotationDelta = currentQuat.clone().multiply(this.prevOrientation.clone().invert());

    // Cache for next frame
    this.prevPosition.copy(currentPos);
    this.prevOrientation.copy(currentQuat);

    // Apply same transform to the locked view
    const camera = this.view.camera;
    const controls = this.view.controls;
    this.applyTransform(translationDelta, rotationDelta, camera, controls);
    controls.update();
  }

  // Applies a rigid-body transform delta (translation + rotation) to the active camera
  // while preserving the user's relative viewpoint (offset from target).
  // This is used to keep the camera locked to a moving/rotating object.
  applyTransform(translation, rotation, camera, controls) {

    // Preserve current offset from target (in world space)
    const targetToCamera = camera.position.clone().sub(controls.target);

    // Rotate the offset vector by the object's rotation delta
    targetToCamera.applyQuaternion(rotation);

    // Move the target by the object's translation delta
    controls.target.add(translation);

    // Reconstruct camera position
    camera.position.copy(controls.target).add(targetToCamera);
  }
};

class View {
  constructor(camera, controls) {
    this.camera = camera;
    this.controls = controls;
    this.initialPosition = camera.position.clone();
    this.initialTarget = controls.target.clone();
    this.lock = null;
  }

  lockTo(marker) {
    // Set a lock on the view so that camera will follow the marker
    this.lock = new Lock(marker, this);
  }

  unlock() {
    this.lock = null;
  }
};

class ViewManager {
  constructor(scene, renderer, startPosition) {
    this.activeIdx = null;
    this.views = [];
    this.scene = scene;
    this.renderer = renderer;

    const mainViewIndex = this.create({ position: startPosition },
                                      { target: new THREE.Vector3(0, 0, 0)});
    this.setActive(mainViewIndex);
    return mainViewIndex;
  }

  get(idx) {
    return this.views[idx];
  }

  getActive() {
    return this.views[this.activeIdx];
  }

  getActiveIndex() {
    return this.activeIdx;
  }

  getAllLocked() {
    return this.views.filter(v => v.lock !== null);
  }

  create(cameraConfig, controlsConfig) {
    const camera = new THREE.PerspectiveCamera(
      cameraConfig.fov || 45,
      innerWidth / innerHeight,
      cameraConfig.near || 0.1,
      cameraConfig.far || 10000
    );

    // orbitControls camera should always be in world coordinates
    this.scene.add(camera);

    const controls = new OrbitControls(camera, this.renderer.domElement);
    controls.target.copy(controlsConfig.target);
    camera.position.copy(cameraConfig.position);
    controls.update(); // Sync controls with initial state
    controls.enabled = false; // New view starts disabled

    // Apply specific properties if they exist in the config
    if (controlsConfig.enablePan !== undefined) controls.enablePan = controlsConfig.enablePan;
    if (controlsConfig.enableZoom !== undefined) controls.enableZoom = controlsConfig.enableZoom;
    if (controlsConfig.maxPolarAngle !== undefined) controls.maxPolarAngle = controlsConfig.maxPolarAngle;

    const view = new View(camera, controls);
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
    console.assert(view.lock === null, "Disposing view with an active lock");
  }

  setActive(viewIndex) {

    const curView = this.getActive();

    // Disable current view
    if (curView)
      curView.controls.enabled = false;

    // Switch active view
    this.activeIdx = viewIndex;

    // Set new view default start position and update controls
    const newView = this.getActive();
    newView.camera.position.copy(newView.initialPosition);
    newView.controls.target.copy(newView.initialTarget);
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
    const newViewIndex = this.create({ position: view.camera.position },
                                     { target: view.controls.target });
    return newViewIndex;
  }

  update() {
    // Update positions of locked cameras
    this.views.forEach(view => {
      if (view.lock)
        view.lock.update();
    });
    // Update controls of active camera
    const activeView = this.getActive();
    activeView.controls.update();
    return activeView.camera;
  }

  // Lock active camera into a geostationary orbit around the selected object
  setOrbit(marker, objectToLock) {

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
    view.lockTo(marker);
  }

  // Setup the camera view for a observer on a planet surface
  createObserver(marker, eyeHeight, lookAhead) {

    // Camera postion and target must be in world coordinates
    const cameraPosWorld = new THREE.Vector3(0, eyeHeight, 0);
    const targetPosWorld = new THREE.Vector3(0, eyeHeight, lookAhead);
    marker.localToWorld(cameraPosWorld);
    marker.localToWorld(targetPosWorld);

    const cameraConfig = {
        position: cameraPosWorld,
        fov: 60,
        //rollAngle: Math.PI / 4
    };
    const controlsConfig = {
        target: targetPosWorld,
        enablePan: false,
        enableZoom: true,
        //maxPolarAngle: Math.PI / 2 + toRadians(2)  // Prevent looking down through the ground
    }

    // Create a new view and set a lock on so that camera will not drift
    const viewIndex = this.create(cameraConfig, controlsConfig);
    return viewIndex;
  }
};

export { ViewManager };
