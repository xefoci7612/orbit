/*
  TODO:
    - When locked reset return to "at lock time" position
    - Look at sign in moon position computation
    - Show sun/raise events with a fading side legend
    - implement view from observer
*/

import * as THREE from 'three';
import { OrbitControls } from 'orbitControls';

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
// CAMERA VIEWS SETUP
// ============================================================================

const ViewManager = {
  activeCamera: null,
  activeControls: null,
  views: [], // An array to hold all our camera/controls pairs
};

// Function to switch views
function setActiveView(viewIndex) {

  if (ViewManager.activeControls)
    ViewManager.activeControls.enabled = false;

  // Set the new active camera and controls
  const newView = ViewManager.views[viewIndex];
  ViewManager.activeCamera = newView.camera;
  ViewManager.activeControls = newView.controls;
  ViewManager.activeControls.enabled = true;

  // Important: Update camera aspect ratio on switch
  ViewManager.activeCamera.aspect = innerWidth / innerHeight;
  ViewManager.activeCamera.updateProjectionMatrix();
}

// Function to create a new view and add it to the manager
function createView(cameraConfig, controlsConfig) {
    const camera = new THREE.PerspectiveCamera(
      cameraConfig.fov || 45,
      innerWidth / innerHeight,
      cameraConfig.near || 0.1,
      cameraConfig.far || 10000
    );
    camera.position.fromArray(cameraConfig.position);

    const controls = new OrbitControls(camera, renderer.domElement);
    controls.minDistance = controlsConfig.minDistance || toUnits(1.2 * EARTH_RADIUS_KM);
    controls.maxDistance = controlsConfig.maxDistance || toUnits(50 * MOON_DISTANCE_KM);
    controls.target.fromArray(controlsConfig.target);
    controls.update(); // Sync controls with initial state

    const view = { camera, controls };
    const viewIndex = ViewManager.views.push(view) - 1;
    return viewIndex;
}

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

// Create the initial camera (this will be view 0)
const cameraDistance = toUnits(2 * MOON_DISTANCE_KM); // From Earth center
const cameraAzimuth = toRadians(180); // Degrees clockwise from +Z axis
const cameraElevation = toRadians(10); // Degrees above horizontal plane
const mainViewIndex = createView({
    position: [
        cameraDistance * Math.cos(cameraElevation) * Math.sin(cameraAzimuth),
        cameraDistance * Math.sin(cameraElevation),
        cameraDistance * Math.cos(cameraElevation) * Math.cos(cameraAzimuth)
    ]
}, {
    target: [0, 0, 0]
});

// Set the initial active view
setActiveView(mainViewIndex);
scene.add(ViewManager.activeCamera);

// Store initial camera position and target
const defaultViewPosition = ViewManager.activeCamera.position.clone();
const defaultViewTarget = ViewManager.activeControls.target.clone();

// Tracks state for an observer on a planet
const Observer = {
    object: null,
    moonVisible: false,
    observerOnEarth: false,
    marker: new THREE.Mesh(
      new THREE.SphereGeometry(toUnits(10)),
      new THREE.MeshBasicMaterial({ color: 0xFFFFFF }) // White
    ),
};

// Tracks state for camera lock-to-target behavior
const Lock = {
  lockedObj: null, // object we currently track
  candidates: [], // objects that can be locked
  prevPosition: new THREE.Vector3(), // obj world position last frame
  prevOrientation: new THREE.Quaternion(), // obj world orientation last frame
  tempVec: new THREE.Vector3(), // helper to avoids re-allocations
  tempQuat: new THREE.Quaternion(), // helper quaternion
  marker: new THREE.Mesh(
    new THREE.SphereGeometry(toUnits(100)), // 100 Km radius
    new THREE.MeshBasicMaterial({ color: 0x00ff00 }) // Lime
  ),
};

// Returns first clickable candidate under mouse, or null
function getClickedCandidate(mouseX, mouseY) {

  const raycaster = new THREE.Raycaster();
  const mouse = new THREE.Vector2();

  // Convert mouse position to normalized device coordinates (-1 to +1)
  mouse.x =  (mouseX / window.innerWidth)  * 2 - 1;
  mouse.y = -(mouseY / window.innerHeight) * 2 + 1;

  // Update the raycaster with the camera and mouse position
  raycaster.setFromCamera(mouse, ViewManager.activeCamera);

  // Find intersections with candidates
  const intersects = raycaster.intersectObjects(Lock.candidates);
  if (intersects.length === 0)
    return null;

  return {
    object: intersects[0].object,
    point: intersects[0].point // point on the surface in world space
  };
}

function unlockCamera() {
  Lock.lockedObj.remove(Lock.marker);
  const mainView = ViewManager.views[mainViewIndex];
  ViewManager.activeControls.minDistance = mainView.controls.minDistance;
  Lock.lockedObj = null;
}

// Orient an observer object on the surface of a parent body.
// Apply the rotation needed to align the observer's local +Y axis
// with the surface normal.
function alignToSurfacePlane(observer) {
  const sourceUp = new THREE.Vector3(0, 1, 0);
  const targetUP = observer.position.clone().normalize(); // position in local coordinates
  observer.quaternion.setFromUnitVectors(sourceUp, targetUP);
}

// Place an observer on an object surface
function setObserver(object, surfacePoint) {

  if (!object.geometry.parameters || object.geometry.parameters.radius === undefined) {
    console.error("Observer must be placed on a THREE.SphereGeometry");
    return;
  }

  Observer.object = object;

  // Show a marker on the clicked point
  const marker = Observer.marker;
  marker.position.copy(surfacePoint); // in world coordinates
  object.worldToLocal(marker.position);
  alignToSurfacePlane(marker); // marker position must be in local coordinates
  object.add(marker); // will move in sync with object

  Observer.observerOnEarth = (object == earth); // before checking for horizon!
  Observer.moonVisible = isAboveHorizon(moonAxis, toUnits(MOON_RADIUS_KM));

  // Set camera on marker along center-surface vector, target at radius X10
  const objCenter = object.getWorldPosition(Lock.tempVec);
  const outwardOffset = surfacePoint.clone().sub(objCenter).multiplyScalar(10);
  const target = surfacePoint.clone().add(outwardOffset);

  // Create a dedicated view for this observer
  const viewIndex = createView({
    position: surfacePoint.toArray()
  }, {
    target: target.toArray(),
    //minDistance: toUnits(1),
  });

  return viewIndex;
}

function removeObserver() {
  if (Observer.object) {
    Observer.object.remove(Observer.marker);
    Observer.object = null;
  }
}

// Locks the camera into a geostationary orbit around the selected object
function lockCameraTo(objectToLock, surfacePoint) {

  if (!objectToLock.geometry.parameters || objectToLock.geometry.parameters.radius === undefined) {
    console.error("Object to lock must be placed on a THREE.SphereGeometry");
    return;
  }

  Lock.lockedObj = objectToLock;
  const camera = ViewManager.activeCamera;
  const controls = ViewManager.activeControls;

  // Initialize the "previous" world state
  Lock.lockedObj.getWorldPosition(Lock.prevPosition);
  Lock.lockedObj.getWorldQuaternion(Lock.prevOrientation);

  // Center camera on target for better UX (user can still pan after)
  const objCenter = Lock.prevPosition.clone();
  controls.target.copy(objCenter);
  const radius = Lock.lockedObj.geometry.parameters.radius
  controls.minDistance = radius * 1.2;

  // Reposition camera along center-surface vector, same distance
  const distanceToCenter = camera.position.distanceTo(objCenter);
  const direction = Lock.tempVec.subVectors(surfacePoint, objCenter).normalize();
  camera.position.copy(direction.multiplyScalar(distanceToCenter)).add(objCenter);

  // Show a marker on the clicked point
  const marker = Lock.marker;
  marker.position.copy(surfacePoint);
  Lock.lockedObj.worldToLocal(marker.position);
  Lock.lockedObj.add(marker);

  controls.update();
}

function updateCameraLock() {

  console.assert(Lock.lockedObj, "Called update lock with no locked object");

  // Get current target transform
  const currentWorldPos = Lock.lockedObj.getWorldPosition(Lock.tempVec);
  const currentWorldQuat = Lock.lockedObj.getWorldQuaternion(Lock.tempQuat);

  // Compute transform delta since last frame
  const translationDelta = currentWorldPos.clone().sub(Lock.prevPosition);
  const rotationDelta = currentWorldQuat.clone().multiply(Lock.prevOrientation.clone().invert());

  // Cache current transform for next frame’s delta
  Lock.prevPosition.copy(currentWorldPos);
  Lock.prevOrientation.copy(currentWorldQuat);

  // Capture the camera's current relationship to its target (this preserves user input)
  const camera = ViewManager.activeCamera;
  const controls = ViewManager.activeControls;
  const targetToCamera = camera.position.clone().sub(controls.target);

  // Apply rotation on the target-to-camera vector
  targetToCamera.applyQuaternion(rotationDelta);

  // Apply translation to target to preserve user pan
  controls.target.add(translationDelta);

  // Recompute camera position from updated target + preserved offset
  camera.position.copy(controls.target).add(targetToCamera);
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
Lock.candidates.push(earth);

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
Lock.candidates.push(moon);

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
  if (!Observer.object)
    return false;

  // Get target's position in the observer frame
  const targetWorldPos = target.getWorldPosition(Lock.tempVec).clone();
  const targetLocalPos = Observer.marker.worldToLocal(targetWorldPos.clone());

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

  if (Lock.lockedObj)
    updateCameraLock();

  // Update the active controls and render with the active camera
  if (ViewManager.activeControls) {
    ViewManager.activeControls.update();
  }

  renderer.render(scene, ViewManager.activeCamera);
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
    this.getClickedCandidate = getClickedCandidate;
    this.lockCameraTo = lockCameraTo;
    this.unlockCamera = unlockCamera;
    this.setObserver = setObserver;
    this.removeObserver = removeObserver;
    this.setActiveView = setActiveView;
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
    setActiveView(mainViewIndex);
    ViewManager.activeCamera.position.copy(defaultViewPosition);
    ViewManager.activeControls.target.copy(defaultViewTarget);
  }
  isLocked(object) {
    return (Lock.lockedObj === object);
  }
  observerActive() {
    return Observer.object !== null;
  }
  cloneView() {
    const currentCam = ViewManager.activeCamera;
    const currentControls = ViewManager.activeControls;
    const newViewIndex = createView({
      position: currentCam.position.toArray()
    }, {
      target: currentControls.target.toArray(),
      minDistance: currentControls.minDistance,
      maxDistance: currentControls.maxDistance
    });
    return newViewIndex;
  }
  resize() {
    const camera = ViewManager.activeCamera;
    camera.aspect = innerWidth / innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(innerWidth, innerHeight);
  }
};

export default Simulation;
