/*
  TODO:
    - When locked reset return to "at lock time" position
    - Look at sign in moon position computation
    - Show sun/raise events with a fading side legend
    - implement view from observer
*/

'use strict';

import * as THREE from 'three';
import { ViewManager } from './view.js';

// ============================================================================
// CONSTANTS AND HELPER FUNCTIONS
// ============================================================================

const EARTH_TEXTURE_URL = 'textures/earth_atmos_2048.jpg';
const MOON_TEXTURE_URL = 'textures/moon_1024.jpg';
const SKY_TEXTURE_URL = [  // In celestial coordinates from https://svs.gsfc.nasa.gov/4851
    'textures/sky_px.png', // Right
    'textures/sky_nx.png', // Left
    'textures/sky_py.png', // Top
    'textures/sky_ny.png', // Bottom
    'textures/sky_pz.png', // Front
    'textures/sky_nz.png'  // Back
];

// Physical properties (in kilometers, converted to scene units)
const EARTH_RADIUS_KM = 6371; // Conventional radius for sphere approx.
const MOON_RADIUS_KM = 1737.53;
const MOON_DISTANCE_KM = 384400;
const EARTH_AXIAL_TILT = 23.44; // degrees
const MOON_AXIAL_TILT = 6.68; // degrees

// Scale conversion constants
const KM_PER_UNIT = EARTH_RADIUS_KM / 10; // Earth radius is 10 units in this scale
const toUnits = km => km / KM_PER_UNIT;
const toRadians = degrees => degrees * Math.PI / 180;

const HORIZON_REFRACTION = toRadians(34.5 / 60); // Standard atmospheric refraction in arcminutes

// Time constants for simulation
const SIM_TIME_SPEED_UP = 60; // 1 simulation hour elapses in 1 minute

const HOUR = 3600; // In SI seconds
const SOLAR_DAY = 24 * HOUR; // Mean Solar Day
const JULIAN_YEAR = 365.25 * SOLAR_DAY;

const SIDERAL_DAY = 23.9344696 * HOUR; // Period of a 360 degrees rotation
const SIDEREAL_YEAR = 365.256363004 * SOLAR_DAY; // For epoch J2000.0

const EARTH_ROTATION_PERIOD = SIDERAL_DAY;
const EARTH_ORBITAL_PERIOD = SIDEREAL_YEAR;
const EARTH_PRECESSION_PERIOD = 25772 * JULIAN_YEAR;
const MOON_ORBITAL_PERIOD = 27.554551 * SOLAR_DAY; // Anomalistic period
const MOON_NODE_PRECESSION_PERIOD = 18.612815932 * JULIAN_YEAR;

// Helpers to calculate date differences
function daysSinceSolstice(targetDate) {
  const year = targetDate.getFullYear();
  // The Summer Solstice in 2025 occurs on June 21 at 03:42 UTC.
  // We use this precise time as our "zero angle" reference point.
  const solstice = new Date(Date.UTC(year, 5, 21, 3, 42, 0));
  return (targetDate - solstice) / (1000 * 60 * 60 * 24); // Convert milliseconds to days
}

// Our reference date is September 10, 2025
const REFERENCE_DATE = new Date('2025-09-10T00:00:00Z');
const DAYS_SINCE_SOLSTICE = daysSinceSolstice(REFERENCE_DATE);

/*
    Init all orbital parameters out of JPL Horizons data for September 10, 2025
    https://ssd.jpl.nasa.gov/horizons/app.html

    Reference frame : Ecliptic of J2000.0

    Origin: The center of the Earth

    Fundamental Plane: The Earth's mean orbital plane (the ecliptic) as it was at the J2000.0 epoch.

    Z-axis: A line perpendicular to this ecliptic plane, pointing north.

    X-axis: This primary direction points towards the vernal equinox, the point where the Sun crosses
            the celestial equator from south to north—as it was at the J2000.0 epoch.

    Y-axis: 90° ahead of the X-axis in the direction of Earth's orbit around the Sun.
*/
const MOON = {
  SemiMajorAxis: toUnits(384821.4581097544), // in km
  EC: 0.05212487669881421,          // eccentricity
  IN: toRadians(5.280240463394533), // inclination relative to ecliptic
  OM: toRadians(347.9713278613414), // longitude of the Ascending Node
   W: toRadians(37.57619864565795), // argument of periapsis (perigee)
  MA: toRadians(352.7811711299613), // mean anomaly
};

// ============================================================================
// SCENE SETUP
// ============================================================================

// Create renderer and scene
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(innerWidth, innerHeight);
document.body.appendChild(renderer.domElement);
const scene = new THREE.Scene();
const cubeTextureLoader = new THREE.CubeTextureLoader();
scene.background = cubeTextureLoader.load(SKY_TEXTURE_URL);

// Ambient light for overall scene illumination
const ambientLight = new THREE.AmbientLight(0xffffff, 0.05); // White light, low intensity
scene.add(ambientLight);

// Directional light to simulate sunlight
const sunLight = new THREE.DirectionalLight(0xffffff, 1);
const sunDistance = toUnits(100 * MOON_DISTANCE_KM);
sunLight.position.set(sunDistance, 0, 0); // on the ecliptic plane +X axis
sunLight.lookAt(0, 0, 0);
scene.add(sunLight);

// Set the default view
const startPosition = (function() {
  const cameraDistance = toUnits(2 * MOON_DISTANCE_KM); // From Earth center
  const cameraAzimuth = toRadians(180); // Degrees clockwise from +Z axis
  const cameraElevation = toRadians(10); // Degrees above horizontal plane
  return new THREE.Vector3(
    cameraDistance * Math.cos(cameraElevation) * Math.sin(cameraAzimuth),
    cameraDistance * Math.sin(cameraElevation),
    cameraDistance * Math.cos(cameraElevation) * Math.cos(cameraAzimuth)
  );
})();
const views = new ViewManager(scene, renderer, startPosition);

// List of objects that can have a locked view on them
const lockables = [];

// Tracks state for an observer on a planet
const Observer = {
    object: null,
    marker:null,
    observerOnEarth: false,
    moonVisible: false,
    viewIndex: null,
    tempVec: new THREE.Vector3(),
};

// Helper for debugging: Add axes helpers to visualize local frames
function addAxesHelper(object, size) {
    const axesHelper = new THREE.AxesHelper(size);
    object.add(axesHelper);
    return axesHelper;
}

// Helper used by observer and orbit locking
function createMarker(radius, color) {
  return new THREE.Mesh(
    new THREE.SphereGeometry(radius),
    new THREE.MeshBasicMaterial({ color })
  );
}

function disposeMarker(marker) {
  marker.geometry.dispose();
  marker.material.dispose();
}

// Add a child object, usually a marker, on a surface of
// a body and align local coordinates of added child
function addAligned(marker, object, surfacePoint) {

  // Surface point must be in world coords
  marker.position.copy(surfacePoint);
  object.worldToLocal(marker.position);
  alignToSurfacePlane(marker);
  object.add(marker);
  marker.updateWorldMatrix(true, false);
}

// Orients the local coordinate frame of an object on the surface of a parent body
function alignToSurfacePlane(marker) {

    // Point on the surface must be in parent coordinates
    const p = marker.position;

    // Calculate Local Y-axis (outward normal from sphere)
    const y_local = p.clone().normalize();

    // Calculate Local X-axis as tangential velocity (z, 0, -x) for Y-axis rotation
    const v = new THREE.Vector3(p.z, 0, -p.x);

    // Detect special case of point exactly on the global Y-axis (pole) and
    // fallback on parent X-axis
    const onY = v.lengthSq() === 0;
    const x_local = onY ? new THREE.Vector3(1, 0, 0) : v.normalize();

    // Local Z-axis is cross product of x_local and y_local for right-handed system
    const z_local = new THREE.Vector3().crossVectors(x_local, y_local);

    // Rotation matrix whose columns are the local basis vectors (x_local, y_local, z_local)
    // This matrix transforms from point's local frame to its parent's (sphere's) frame.
    const rotationMatrix = new THREE.Matrix4();
    rotationMatrix.makeBasis(x_local, y_local, z_local);

    // Set point's quaternion from this rotation matrix
    marker.quaternion.setFromRotationMatrix(rotationMatrix);
}

// Returns first clickable candidate under mouse, or null
function getClickedCandidate(mouseX, mouseY) {

  const raycaster = new THREE.Raycaster();
  const mouse = new THREE.Vector2();

  // Convert mouse position to normalized device coordinates (-1 to +1)
  mouse.x =  (mouseX / window.innerWidth)  * 2 - 1;
  mouse.y = -(mouseY / window.innerHeight) * 2 + 1;

  // Update the raycaster with the camera and mouse position
  const view = views.getActive();
  raycaster.setFromCamera(mouse, view.camera);

  // Find intersections with candidates
  const intersects = raycaster.intersectObjects(lockables);
  if (intersects.length === 0)
    return null;

  return {
    object: intersects[0].object,
    point: intersects[0].point // point on the surface in world space
  };
}

// Locks the camera into a geostationary orbit around the selected object
function lockCameraTo(objectToLock, surfacePoint) {

  // Set a marker in local coordinates on object's surface point
  const marker = createMarker(toUnits(50), 0x00ff00);
  addAligned(marker, objectToLock, surfacePoint);

  // Lock the active view aligned to object's center
  views.setOrbit(marker, objectToLock);
}

// Place an observer on an object surface
function setObserver(object, surfacePoint) {

  // Currently we handle only one observer
  console.assert(Observer.object === null, "Setting already existing observer");

  // Position and orient the marker on the surface
  const marker = createMarker(toUnits(10), 0xFFFFFF);
  addAligned(marker, object, surfacePoint);

  Observer.object = object;
  Observer.marker = marker;

  // Init observer events
  Observer.observerOnEarth = (object === earth);
  Observer.moonVisible = isAboveHorizon(moonAxis, toUnits(MOON_RADIUS_KM));

  // Create a new observer view placed on marker
  const eyeHeight = toUnits( 20); // km above surface
  const lookAhead = toUnits(100); // look at km ahead on horizon
  const viewIndex = views.createObserver(marker, eyeHeight, lookAhead);
  Observer.viewIndex = viewIndex;

  // Lock the view camera on the marker
  const view = views.get(viewIndex);
  view.lockTo(marker);

  // DEBUG
  addAxesHelper(marker, 10);
  addAxesHelper(view.camera, 5);

  return viewIndex;
}

function removeObserver() {
  const view = views.get(Observer.viewIndex);
  view.unlock();
  const marker = Observer.marker;
  marker.parent.remove(marker);
  disposeMarker(marker);
  views.dispose(Observer.viewIndex);
  Observer.object = null;
  Observer.viewIndex = null;
  Observer.marker = null;
}

// ============================================================================
// CELESTIAL BODY CREATION
// ============================================================================

/* WARNING: By default Three.js applies rotations in XYZ order. We use nested
   Object3D containers to apply rotations in a different order.

   Three.js coordinate system is +X right, +Y up, +Z toward camera
*/

// Computes an array of 3D vertex on an elliptic curve in XZ plane with
// Earth at (0,0,0). The local X-axis represents the "line of nodes,"
// and the positive direction (+X, 0, 0) points directly at the ascending node.
function ellipticCurve(SemiMajorAxis, EC, W) {
  const SEGMENTS = 360;
  const points = new Float32Array(SEGMENTS * 3);

  // Distance from Earth given true anomaly ta (ellipse polar equation)
  const distance = ta => SemiMajorAxis * (1 - EC * EC) / (1 + EC * Math.cos(ta));

  // Loop starts at perigee (W) when i = 0, and angle 0 (ta = -W)
  // corresponds to the ascending node.
  for (let i = 0; i < SEGMENTS; i++) {
    const ta = (i / SEGMENTS) * 2 * Math.PI; // true anomaly
    const r = distance(ta); // distance from Earth at true anomaly
    const angle = ta + W; // angle from ascending node
    points[i*3 + 0] = r * Math.cos(angle); // x
    points[i*3 + 1] = 0;                   // y  (on orbital plane)
    points[i*3 + 2] = r * Math.sin(angle); // z
  }

  // Find ascending and descending nodes x-axis coordinates
  const rAsc = distance(-W);
  const rDes = distance(-W + Math.PI);
  return [points, rAsc, rDes];
}

// Create Earth hierarchy
// The innermost node performs daily rotations around pole axis
// Earth texture map is already centered on 0° longitude.
const textureLoader = new THREE.TextureLoader();
const earthTexture = textureLoader.load(EARTH_TEXTURE_URL);
const earth = new THREE.Mesh(
  new THREE.SphereGeometry(toUnits(EARTH_RADIUS_KM), 64, 64),
  new THREE.MeshStandardMaterial({ map: earthTexture, roughness: 0.8 })
);
lockables.push(earth);

// Add parent to tilt and rotate vertical axis for precession
const earthAxis = new THREE.Object3D();
// Set the Earth's obliquity (axial tilt).
// At summer solstice the North Pole must be tilted towards the Sun.
// In Three.js coordinate system (+X right, +Y up, +Z toward camera),
// when Sun is on the +X axis, a negative (clockwise) rotation around
// the Z-axis correctly tilts the Earth's North Pole towards the sun.
earthAxis.rotation.z = -toRadians(EARTH_AXIAL_TILT);
earthAxis.add(earth);

// Create Moon hierarchy, each node performs at most one rotation
// Starting from the innermost node, the moon mesh that performs
// rotation around its axis. Texture map is already centered on
// 0° longitude.
const moonTexture = textureLoader.load(MOON_TEXTURE_URL);
const moon = new THREE.Mesh(
  new THREE.SphereGeometry(toUnits(MOON_RADIUS_KM), 32, 32),
  new THREE.MeshStandardMaterial({ map: moonTexture, roughness: 0.9 })
);

// Add parent to tilt and rotate vertical axis for precession
const moonAxis = new THREE.Object3D();
moonAxis.rotation.z = toRadians(MOON_AXIAL_TILT); // FIXME initial offset missing
moonAxis.add(moon);
lockables.push(moon);

// Add parent to set a fixed inclined orbital path, this is the target when setting
// current moon position in animation (through moonAxis object)
const moonOrbit = new THREE.Object3D();
moonOrbit.rotation.x = MOON.IN; // inclination to ecliptic
moonOrbit.add(moonAxis);

// Create an elliptical curve for the Moon's orbit in its orbital plane
const [points, rAsc, rDes] = ellipticCurve(MOON.SemiMajorAxis, MOON.EC, MOON.W);
const vertexBuffer = new THREE.BufferAttribute(points, 3);
const orbitGeometry = new THREE.BufferGeometry().setAttribute('position', vertexBuffer);
const orbitMaterial = new THREE.LineBasicMaterial({
  color: 0x4488ff, // light blue
  transparent: true,
  opacity: 0.35
});
const orbitPath = new THREE.LineLoop(orbitGeometry, orbitMaterial);
moonOrbit.add(orbitPath);

// Setup the node markers corresponding to ascending/descending nodes
function getMarker(color, position) {
  const geometry = new THREE.SphereGeometry(toUnits(600));
  const material = new THREE.MeshBasicMaterial({ color: color });
  const marker = new THREE.Mesh(geometry, material);
  marker.position.set(...position);
  return marker;
}
// Nodes lay on x-axis in local coordinates
moonOrbit.add(getMarker(0xff00ff, [+rAsc, 0, 0])); // Magenta
moonOrbit.add(getMarker(0xffff00, [-rDes, 0, 0])); // Yellow

// Add parent to follow the Moon's orbit nodal precession (swivel)
const ascendingNodePivot = new THREE.Object3D();
ascendingNodePivot.add(moonOrbit);

// Finally add Moon and Earth hierarchy to scene
scene.add(ascendingNodePivot);
scene.add(earthAxis);

// ============================================================================
// CELESTIAL BODY RUNTIME COMPUTATIONS
// ============================================================================

/**
 * Solves Kepler's Equation M = E - e * sin(E) for the Eccentric Anomaly (E).
 * Uses the Newton-Raphson iterative method.
 *
 * @param {number} M - The Mean Anomaly in radians.
 * @param {number} e - The eccentricity of the orbit.
 * @returns {number} The Eccentric Anomaly (E) in radians.
 */
function solveKepler(M, e) {
  let E = M; // Initial guess for E is M
  const maxIterations = 10;
  const tolerance = 1e-6; // A small tolerance for accuracy

  for (let i = 0; i < maxIterations; i++) {
    const f = E - e * Math.sin(E) - M; // The function we want to find the root of
    const fPrime = 1 - e * Math.cos(E); // The derivative df/dE
    const deltaE = f / fPrime;

    E = E - deltaE; // New, better guess for E

    if (Math.abs(deltaE) < tolerance)
      break; // Solution is accurate enough
  }
  return E;
}

function getMoonPosition(elapsedTime) {

  // Calculate current Mean Anomaly (M) and Eccentric Anomaly (E)
  const n = 2 * Math.PI / MOON_ORBITAL_PERIOD; // Mean motion (rad/sec)
  const M = MOON.MA - n * elapsedTime; // Current Mean Anomaly
  const EC = MOON.EC;
  const E = solveKepler(M, EC); // Kepler's Second Law

  // Calculate the current True Anomaly (TA or ν), the real angle from perigee
  const [sinE, cosE] = [Math.sin(E), Math.cos(E)];
  const sinTA = Math.sqrt(1 - EC*EC) * sinE / (1 - EC * cosE);
  const cosTA = (cosE - EC) / (1 - EC * cosE);
  const TA = Math.atan2(sinTA, cosTA);

  // Calculate distance from Earth (focus)
  const r = MOON.SemiMajorAxis * (1 - EC * EC) / (1 + EC * Math.cos(TA));

  // Calculate position in the orbital plane
  const angleFromAscendingNode = TA + MOON.W; // x-axis is at ascending node
  const x = r * Math.cos(angleFromAscendingNode);
  const z = r * Math.sin(angleFromAscendingNode);
  return [x, 0, z]; // Return 3D coordinates (Y=0 in orbital plane)
}

function isAboveHorizon(target, radius) {

  // Get target's position in the observer frame
  const temp = Observer.tempVec;
  target.getWorldPosition(temp);
  const targetLocalPos = Observer.marker.worldToLocal(temp);

  // Correction for the object's upper limb (its radius)
  const radiusCorrection = radius;

  // Correction for atmospheric refraction (lifts the image)
  let refractionCorrection = 0;
  if (Observer.observerOnEarth) {
    const distance = targetLocalPos.length();
    refractionCorrection = distance * Math.tan(HORIZON_REFRACTION);
  }

  // The object is visible if its center is above this negative threshold
  const visibilityThreshold = -(radiusCorrection + refractionCorrection);
  return targetLocalPos.y > visibilityThreshold;
}

// ============================================================================
// ANIMATION LOOP
// ============================================================================

// Calculate simulation constants and state
const earthRotationSpeed = 2 * Math.PI / EARTH_ROTATION_PERIOD;
const precessionSpeed = 2 * Math.PI / EARTH_PRECESSION_PERIOD;
const moonSiderealRotationSpeed = 2 * Math.PI / MOON_ORBITAL_PERIOD;
const moonNodeSpeed = 2 * Math.PI / MOON_NODE_PRECESSION_PERIOD;
const sunOrbitSpeed = 2 * Math.PI / EARTH_ORBITAL_PERIOD;
const sunInitialAngle = (DAYS_SINCE_SOLSTICE / 365.25) * 2 * Math.PI;
// Prime meridian is at midnight at solstice 00:00
const earthInitialAngle = sunInitialAngle + Math.PI;

// Clock works also when tab is hidden, can be set/reset by UI, time
// can run faster and/or backward from real time and is adjustable
const Clock = class {
  constructor() {
    this.running = true;
    this.referenceTime = REFERENCE_DATE.getTime(); // fixed t=0 point
    this.reset();
  }
  reset() {
    this.simulationTime = this.referenceTime;
    this.realTimeAtResume = Date.now();
    this.speedX = SIM_TIME_SPEED_UP;
  }
  _deltaSinceLastResume() {
    // Accumulate time only while running. It is simulation time, not real time
    return this.running ? this.speedX * (Date.now() - this.realTimeAtResume) : 0;
  }
  elapsedTime() {
    const simTime = this.simulationTime + this._deltaSinceLastResume();
    return (simTime - this.referenceTime) / 1000; // in secs
  }
  togglePause() {
    if (this.running) {
      // Before stopping, add the "ticking" time to the "banked" time
      this.simulationTime += this._deltaSinceLastResume();
    } else {
      this.realTimeAtResume = Date.now();
    }
    this.running = !this.running;
    return this.running;
  }
  speed() {
    return this.speedX / SIM_TIME_SPEED_UP; // SIM_TIME_SPEED_UP -> 1
  }
  setSpeed(speed) {
    // Banks any accumulated time at the old speed by simulating a
    // stop/start, then sets the new speed. It works whether the
    // clock is running or paused.
    this.togglePause();
    this.togglePause();
    this.speedX = SIM_TIME_SPEED_UP * speed; // 1 -> SIM_TIME_SPEED_UP
  }
  setDate(newDate) {
    // Adjust the "banked" time to make the total time equal the newDate
    this.simulationTime = newDate.getTime() - this._deltaSinceLastResume();
  }
};
const simClock = new Clock();

const simEvents = { atRise: false, atSet: false };

// Animation loop function
function animate() {

  const elapsedTime = simClock.elapsedTime(); // in seconds

  // Rotate Earth
  earth.rotation.y = earthInitialAngle  + elapsedTime * earthRotationSpeed;

  // Earth axial precession (at t = 0 inclination is aligned toward Sun at solstice)
  earthAxis.rotation.y = precessionSpeed * elapsedTime;

  // Moon rotation around its axis, in tidal locking with Earth
  // Offset 180-degree because texture is centered on 0° longitude (facing earth)
  moon.rotation.y = Math.PI + moonSiderealRotationSpeed * elapsedTime;

  // Moon orbit around Earth
  const [x, y, z] = getMoonPosition(elapsedTime);
  moonAxis.position.set(x, y, z);

  // The Moon's orbital plane precesses backwards (retrograde).
  // This is a rotation around the Ecliptic Pole (the Y-axis).
  ascendingNodePivot.rotation.y = MOON.OM - moonNodeSpeed * elapsedTime;

  // The sun's orbital motion is counter-clockwise when viewed from above (positive Y)
  // In Three.js coordinate system (+X right, +Y up, +Z toward camera).
  // So a counter-clockwise rotation requires the sun to move from positive +X toward
  // negative -Z and angle will decrease with time.
  const sunAngle = -(sunInitialAngle + elapsedTime * sunOrbitSpeed);
  sunLight.position.set(sunDistance * Math.cos(sunAngle), 0, sunDistance * Math.sin(sunAngle));

  if (Observer.object) {
    // Pause animation at Moon raise/set
    const isMoonVisible = isAboveHorizon(moonAxis, toUnits(MOON_RADIUS_KM));
    const atRaise =  isMoonVisible && !Observer.moonVisible;
    const atSet   = !isMoonVisible &&  Observer.moonVisible;
    Observer.moonVisible = isMoonVisible;
    Object.assign(simEvents, { atRise: atRaise, atSet: atSet });
  }

  // Update the active controls and render with the active camera
  const activeCamera = views.update();
  renderer.render(scene, activeCamera);
  return simEvents;
}

// ============================================================================
// EXPORT CLASS
// ============================================================================

class Simulation {
  constructor() {
    this.update = animate;
    this.setDate = simClock.setDate.bind(simClock);
    this.speed = simClock.speed.bind(simClock);
    this.setSpeed = simClock.setSpeed.bind(simClock);
    this.togglePause = simClock.togglePause.bind(simClock);
    this.setActiveView = views.setActive.bind(views);
    this.disposeView = views.dispose.bind(views);
    this.lockCameraTo = lockCameraTo;
    this.setObserver = setObserver;
    this.removeObserver = removeObserver;
    this.getClickedCandidate = getClickedCandidate;
  }
  getRenderer() {
    return renderer;
  }
  getDate() {
    const elapsedMsec = simClock.elapsedTime() * 1000;
    const date = new Date(REFERENCE_DATE.getTime() + elapsedMsec);
    return [date, elapsedMsec];
  }
  reverseSpeed() {
    const speed = -simClock.speed();
    simClock.setSpeed(speed);
    return speed;
  }
  reset() {
    simClock.reset();
    views.setDefault();
  }
  isLocked(object) {
    const lockedViews = views.getAllLocked();
    return lockedViews.some(v => v.lock.object === object);
  }
  unlockCamera() {
    const view = views.getActive();
    const marker = view.lock.target;
    marker.parent.remove(marker);
    disposeMarker(marker);
    view.unlock();
  }
  cloneView() {
    const v = views.getActive();
    const newViewIndex = views.clone(v);
    return newViewIndex;
  }
  resize() {
    const v = views.getActive();
    const camera = v.camera;
    camera.aspect = innerWidth / innerHeight;
    camera.updateProjectionMatrix();
    views.renderer.setSize(innerWidth, innerHeight);
  }
};

export default Simulation;
