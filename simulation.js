'use strict';

/*
  TODO:
    - When locked reset return to "at lock time" position
    - Show sun/raise events with a fading side legend
*/


/*
  Our World Frame is a right-handed coordinate system derived from ICRS Heliocentric
  Ecliptic J2000 Frame and adapted to Three.js conventions for axes

  ICRS Heliocentric Ecliptic J2000 Frame

    Ecliptic J2000 Frame refers to Heliocentric Ecliptic coordinates according to the
    International Celestial Reference System (ICRS) at the J2000.0 epoch.
    It is a heliocentric, right-handed coordinate system.

    EPOCH (J2000.0):
        The reference frame is fixed to a specific moment: January 1, 2000, at 12:00
        Terrestrial Time (TT), which is Julian Date 2451545.0 TT.

    ORIGIN:
        The formal origin (0,0,0) is the Barycenter of the Solar System (SSB)

    THE FUNDAMENTAL PLANE (THE XY PLANE):
        This is the Earth's mean orbital plane as it existed at the J2000.0 epoch,
        also known as the "Ecliptic of J2000.0"

    AXES DEFINITION:
     +X Axis: The principal direction. It points from the Solar System Barycenter
              towards the mean Vernal Point as defined at the J2000.0 epoch.
              The Vernal Point (First Point of Aries) is one of the two points at the
              intersection of the celestial equator and the ecliptic, and its direction
              is defined by the Earth-to-Sun vector at the instant of the Vernal Equinox.

      +Z Axis: Points towards the North Ecliptic Pole at the J2000.0 epoch. It is
               perpendicular to the XY plane

      +Y Axis: Completes the right-handed coordinate system (X × Y = Z). This axis
               points towards the Winter Solstice direction. At the Vernal Equinox,
               the +Y axis has opposite direction of Earth's velocity.


  Scene World Coordinates Frame

  We use a right-handed heliocentric coordinate system (referred as World Coordinates)
  adapted to Three.js conventions for axes.

  Differences from ICRS Heliocentric Ecliptic J2000 Frame:

     ORIGIN: The formal origin (0,0,0) is the center of the Sun

    +Y Axis: Is the +Z Axis in ICRS Frame (UP direction)

    +Z Axis: Is the -Y Axis in the ICRS Frame. It points toward Summer Solstice, at
             Vernal Equinox time, the +Z axis has same direction of Earth's velocity.
*/

import * as THREE from 'three';
import { ViewManager } from './view.js';
import { CelestialEventManager } from './events.js';
import { Observer } from './observer.js';
import { init_debug, set_sunline_length, addAxesHelper } from './validate.js';

export const DEBUG = false; // Enable axes and objects visualizations for debug purposes
export const VALIDATE = false; // Dry-run validation instead of normal visualization mode

const EARTH_TEXTURE_URL = 'textures/earth_atmos_2048.jpg';
const MOON_TEXTURE_URL = 'textures/moon_1024.jpg';

// In celestial coordinates from https://svs.gsfc.nasa.gov/4851
// FIXME convert to Heliocentric Ecliptic ICRS Frame
const SKY_TEXTURE_URL = [
    'textures/sky_px.png', // Right
    'textures/sky_nx.png', // Left
    'textures/sky_py.png', // Top
    'textures/sky_ny.png', // Bottom
    'textures/sky_pz.png', // Front
    'textures/sky_nz.png'  // Back
];

// Physical conventional properties in kilometers or degrees
const EARTH_RADIUS_KM  = 6371;
const MOON_RADIUS_KM   = 1737.53;
const SUN_RADIUS_KM    = 695700;
const MOON_DISTANCE_KM = 384400;
const EARTH_AXIAL_TILT = 23.439291; // degrees
const MOON_AXIAL_TILT  = 6.68;      // degrees

// Earth-Moon Barycenter (EBM) Correction Constant to convert
// from EMB->Moon distance to EMB->Earth distance
// Barycenter rule is:
//    m₁ × d₁ = m₂ × d₂  --> emb_e = emb_m * MOON_EARTH_MASS_RATIO
const MOON_EARTH_MASS_RATIO = 1 / 81.3;

// The official IAU 2012 definition of the Astronomical Unit
const AU = 149597870.7; // in kilometers

// Scale conversions
const KM_PER_UNIT = EARTH_RADIUS_KM / 50; // Earth radius to units
const toUnits = km => km / KM_PER_UNIT;
const fromUnits = u => u * KM_PER_UNIT;
const toRadians = degrees => degrees * Math.PI / 180;
const toDegrees = rad => rad * 180 / Math.PI;
const HMSToRadians = s => s.split(/\s+/).reduce((acc, v, i) => acc + parseFloat(v) / [1, 60, 3600][i], 0) * (Math.PI / 12);

// Helper to normalize an angle in radians to the [0, 2PI] range
const normalizeRad = rad => {
  let normalized = rad % (2 * Math.PI);
  if (normalized < 0) { normalized += 2 * Math.PI; }
  return normalized;
};

// Rescaled Sun distance for placing directional Sun light source
const SUN_LIGHT_DISTANCE = toUnits(10 * MOON_DISTANCE_KM);

// Time constants for simulation
const SIM_TIME_SPEED_UP = 60; // 1 simulation hour elapses in 1 minute

// Reference Epoch for simulation is J2000.0
const J2000_EPOCH = new Date('2000-01-01T12:00:00Z');

const HOUR = 3600; // In SI seconds
const SOLAR_DAY = 24 * HOUR; // Mean Solar Day
const SIDERAL_DAY = 23.9344696 * HOUR; // Period of a 360 degrees rotation
const JULIAN_YEAR = 365.25 * SOLAR_DAY;

const EARTH_AXIAL_PERIOD = 25772 * JULIAN_YEAR;
const MOON_ORBITAL_PERIOD = 27.32166 * SOLAR_DAY;

const GAST = HMSToRadians("08 52 45.4526"); // Greenwich Apparent Sidereal Time (GAST)
const GAST_EPOCH = new Date('2025-09-10T09:34:03Z'); // Time of when GAST was sampled
const GAST_TIME = (GAST_EPOCH.getTime() - J2000_EPOCH.getTime()) / 1000; // in secs

// Standard atmospheric refraction in arcminutes
const HORIZON_REFRACTION = toRadians(34.5 / 60);

/*
    Earth Mean Orbital elements

    Reference plane: Ecliptic plane at J2000 Epoch

    Source: https://www.met.reading.ac.uk/~ross/Astronomy/Planets.html
*/

// Target Body: Earth, Coordinate Center: Sun
const EARTH_J2000 = {
  A: toUnits(1.00000011 * AU), // Semi-major axis in AU
  EC: 0.01671022,              // Eccentricity
  IN: toRadians(0.00005),      // Inclination to Ecliptic
  // For near-zero inclinations, Longitude of Ascending Node (OM) is
  // ill-defined and Longitude of Perihelion (W_bar or ϖ) is used
  // instead of Argument of Periapsis.
  W_bar: toRadians(102.94719), // Longitude of Perihelion (ϖ)
  L: toRadians(100.46435),     // Mean Longitude (L = ϖ + M)
};

// Rate of Change (per Julian century)
const EARTH_CHANGE_RATES = {
  A: toUnits(-0.00000005 * AU),
  EC: -0.00003804,
  IN: toRadians(-46.94 / 3600), // From arcsecs, change relative to the ecliptic of the date
  W_bar: toRadians(1198.28 / 3600), // From arcsecs, rate for Longitude of Perihelion (ϖ)
  L: toRadians(129597740.63 / 3600), // From arcsecs, rate for Mean Longitude (L)
}

// Propagate Earth parameters since reference epoch. For earth update is
// based on a fixed and very simple set of change rates.
function propagateEarth(julianCenturies) {
  const current = {};
  for (const key of Object.keys(EARTH_J2000)) {
    const baseValue = EARTH_J2000[key];
    current[key] = baseValue + EARTH_CHANGE_RATES[key] * julianCenturies;

    // The mean anomaly and other parameters can grow to a very large number.
    // It's good practice to normalize it to the 0-2PI range.
    if (key !== 'A' && key !== 'EC') {
      current[key] = normalizeRad(current[key]);
    }
  };

  return current;
}

/*
  The following polynomial expressions for the Moon's mean elements are based on
  a simplified version of the ELP-2000/82B lunar theory.

  Source: "Astronomical Algorithms" (2nd Edition) by Jean Meeus, Chapter 47.
*/
function propagateMoon(julianCenturies) {

  const T = julianCenturies;

  // Moon's Mean Longitude (L') in degrees
  const L_prime = 218.3164477 + 481267.88123421 * T - 0.0015786 * T*T + T*T*T/538841 - T*T*T*T/65194000;

  // Moon's Mean Anomaly (M')
  const M_prime = 134.9633964 + 477198.8675055 * T + 0.0087414 * T*T + T*T*T/69699 - T*T*T*T/14712000;

  // Moon's Argument of Latitude (F)
  const F = 93.2720950 + 483202.0175233 * T - 0.0036539 * T*T - T*T*T/3526000 + T*T*T*T/863310000;

  // Derive the required orbital elements from the fundamental arguments
  const meanAnomalyDeg = M_prime;
  const argumentOfPeriapsisDeg = F - M_prime;
  const longAscendingNodeDeg = L_prime - F;

  return {
    A:  toUnits(379700),
    EC: 0.0554,
    IN: toRadians(5.145),
    MA: normalizeRad(toRadians(meanAnomalyDeg)),
    W:  normalizeRad(toRadians(argumentOfPeriapsisDeg)),
    OM: normalizeRad(toRadians(longAscendingNodeDeg)),
  };
}

// Example: ISS Orbital Elements (Osculating, for a specific epoch)
// Epoch: 2023-10-27 00:00:00 UTC
const ISS_EPOCH = new Date('2023-10-27T00:00:00Z');
const ISS_TIME = (ISS_EPOCH.getTime() - J2000_EPOCH.getTime()) / 1000; // in secs

// Convert km to scene units
const ISS_J2000 = {
  A:  toUnits(EARTH_RADIUS_KM + 417), // Semi-major axis (km) -> units
  EC: 0.0007,                        // Eccentricity
  IN: toRadians(51.64),              // Inclination (deg) -> rad
  MA: toRadians(130.5),              // Mean Anomaly (deg) -> rad
  W:  toRadians(247.4),              // Argument of Perigee (deg) -> rad
  OM: toRadians(325.0),              // RAAN (deg) -> rad
  EPOCH: ISS_EPOCH
};

// Earth-specific constants for J₂ propagation
const EARTH_MU = 398600.4418; // Standard gravitational parameter (km³/s²)
const EARTH_J2 = 1.08262668e-3; // J₂ coefficient (dimensionless)

// Use a J₂ Perturbation Propagator to update satellite position
function propagateSatellite(initialElements, timeSinceEpoch) {

  // a, e, and i are constant in this model
  const current = {
     A: initialElements.A,
    EC: initialElements.EC,
    IN: initialElements.IN,
  };

  const a  = fromUnits(initialElements.A); // Semi-major axis in km
  const e  = initialElements.EC;           // Eccentricity
  const i  = initialElements.IN;           // Inclination in radians
  const M0 = initialElements.MA;           // Mean anomaly at epoch
  const w0 = initialElements.W;            // Arg of Perigee at epoch
  const O0 = initialElements.OM;           // RAAN at epoch

  // Mean motion (n) in rad/sec
  const n = Math.sqrt(EARTH_MU / Math.pow(a, 3));

  // Semi-latus rectum (p) in km
  const p = a * (1 - e * e);

  // Pre-compute the J₂ perturbation factor for efficiency
  const j2_factor = -(3 / 2) * n * EARTH_J2 * Math.pow(EARTH_RADIUS_KM / p, 2);

  // Rate of change for RAAN (Ω̇)
  const raan_dot = j2_factor * Math.cos(i);

  // Rate of change for Argument of Perigee (ω̇)
  const argp_dot = j2_factor * (2 - (5 / 2) * Math.sin(i) * Math.sin(i));

  // Update RAAN
  current.OM = normalizeRad(O0 + raan_dot * timeSinceEpoch);

  // Update Argument of Perigee
  current.W = normalizeRad(w0 + argp_dot * timeSinceEpoch);

  // Update Mean Anomaly
  current.MA = normalizeRad(M0 + n * timeSinceEpoch);

  return current;
}

// ============================================================================
// SCENE SETUP
// ============================================================================

/*
  Setup a hierarchy of 3D Objects to split complex movement in independent
  rotations
*/

const temp1Vec = new THREE.Vector3();
const temp2Vec = new THREE.Vector3();
const temp3Vec = new THREE.Vector3();
const earthPosWorldVec = new THREE.Vector3();
const moonPosWorldVec = new THREE.Vector3();
const embPosWorldVec = new THREE.Vector3();
const toWorldQuat = new THREE.Quaternion();

// Create a marker object
function getMarker(radius, color, position) {
  const geometry = new THREE.SphereGeometry(radius);
  const material = new THREE.MeshBasicMaterial({ color: color });
  const marker = new THREE.Mesh(geometry, material);
  marker.position.copy(position);
  return marker;
}

// Release marker resources
function disposeMarker(marker) {
  marker.geometry.dispose();
  marker.material.dispose();
}

// Create an ellipses path line, local +X axis points toward periapsis
function getOrbitPath(A, EC) {
  const points = ellipticCurve(A, EC);
  const vertexBuffer = new THREE.BufferAttribute(points, 3);
  const orbitPathGeo = new THREE.BufferGeometry().setAttribute('position', vertexBuffer);
  const orbitPathMat = new THREE.LineBasicMaterial({
    color: 0x4488ff, // light blue
    transparent: true,
    opacity: 0.35
  });
  const orbitPath = new THREE.LineLoop(orbitPathGeo, orbitPathMat);
  return orbitPath;
}

// Create renderer and scene
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(innerWidth, innerHeight); // Full viewport
document.body.appendChild(renderer.domElement);
const textureLoader = new THREE.TextureLoader();
const cubeTextureLoader = new THREE.CubeTextureLoader();
const scene = new THREE.Scene();
scene.background = cubeTextureLoader.load(SKY_TEXTURE_URL);

// Ambient light for overall scene illumination
const ambientLight = new THREE.AmbientLight(0xffffff, 0.05); // White light, low intensity
scene.add(ambientLight);

// Directional light for sun light, set at the Origin
const sunLight = new THREE.DirectionalLight(0xffffff, 1);
sunLight.position.set(0, 0, 0); // default is 'looking from top'
scene.add(sunLight);
scene.add(sunLight.target); // target must be added too!


// Sun hierarchy
//
// The Earth/EMB 2D (apsidal) orbital plane
// Sun is at the origin, +X axis points toward Perihelion
const earthOrbitalPlane = new THREE.Object3D();
scene.add(earthOrbitalPlane);

// EMB (Earth-Moon Barycenter) pivot. It is positioned within the
// earthOrbitalPlane to follow the elliptical path. An inverse
// rotation is applied to it, cancelling the parent's rotation.
// The result is that this pivot's local axes remain aligned with
// the World axes, providing a stable reference for Earth's axial tilt.
const embPivot = new THREE.Object3D();
earthOrbitalPlane.add(embPivot);


// Moon hierarchy
//
// The Moon nodes orbital plane, it is set as child of scene because
// orbital parameters (like inclinations) refer to J2000 Epoch
// In animation we transalte the plane to set EMB at the origin and
// rotate Y axis so that local +X axis points toward Ascending Node
const moonNodesOrbitalPlane = new THREE.Object3D();
scene.add(moonNodesOrbitalPlane);

// Add markers on the Moon's orbit ascending and descending nodes
const ascendingNodeMarker  = getMarker(toUnits(600), 0xff00ff, temp1Vec); // Magenta
const descendingNodeMarker = getMarker(toUnits(600), 0xff00ff, temp1Vec); // Yellow
moonNodesOrbitalPlane.add(ascendingNodeMarker);
moonNodesOrbitalPlane.add(descendingNodeMarker);

// The Moon absidal orbital plane
// EMB is at the origin, same XZ plane as the nodes orbital plane,
// local +X axis points toward periapsis / perigee
const moonAbsidalOrbitalPlane = new THREE.Object3D();
moonNodesOrbitalPlane.add(moonAbsidalOrbitalPlane);

// Moon with tilted vertical axis, it is the object that orbits
// around EMB, local +X axis points toward periapsis / perigee
const tiltedMoon = new THREE.Object3D();
moonAbsidalOrbitalPlane.add(tiltedMoon);

// The actual Moon Mesh object, in tidal locked rotation
// around its pole axis
const moonTexture = textureLoader.load(MOON_TEXTURE_URL);
const moon = new THREE.Mesh(
  new THREE.SphereGeometry(toUnits(MOON_RADIUS_KM), 32, 32),
  new THREE.MeshStandardMaterial({ map: moonTexture, roughness: 0.9 })
);
tiltedMoon.add(moon);

// An elliptical curve to show Moon's orbit, local +X axis points
// toward periapsis / perigee
// Semi axis major and eccentricity are assumed to be constant
const moonEllipses = { A: toUnits(379700), EC: 0.0554 };
const moonOrbitPath = getOrbitPath(moonEllipses.A, moonEllipses.EC);
moonAbsidalOrbitalPlane.add(moonOrbitPath);


// Earth hierarchy
//
// Earth with tilted vertical axis, it is the object that orbits
// around EMB point.
// In animation is set at a computed offset from EMB, local +X axis
// is always aligned with World +X axis (toward Vernal Point)
const tiltedEarth = new THREE.Object3D();
embPivot.add(tiltedEarth);

// The actual Earth Mesh object, performs daily rotation
// around its pole axis
const earthTexture = textureLoader.load(EARTH_TEXTURE_URL);
const earth = new THREE.Mesh(
  new THREE.SphereGeometry(toUnits(EARTH_RADIUS_KM), 64, 64),
  new THREE.MeshStandardMaterial({ map: earthTexture, roughness: 0.8 })
);
tiltedEarth.add(earth);

// Satellite hierarchy
//
// The entire satellite system is parented to the Earth
const satelliteNodesOrbitalPlane = new THREE.Object3D();
tiltedEarth.add(satelliteNodesOrbitalPlane);

const satelliteAbsidalOrbitalPlane = new THREE.Object3D();
satelliteNodesOrbitalPlane.add(satelliteAbsidalOrbitalPlane);

// The actual Satellite Mesh object
const SAT_VISIBILITY_SCALE = 20; // in km
const SAT_REALISTIC_SCALE = 2 / SAT_VISIBILITY_SCALE;
const satellite = new THREE.Mesh(
  new THREE.SphereGeometry(toUnits(SAT_VISIBILITY_SCALE), 8, 8),
  new THREE.MeshBasicMaterial({ color: 0xff0000 }) // Red
);
satelliteAbsidalOrbitalPlane.add(satellite);

// An elliptical curve to show satellite's orbit, local +X axis points
// toward periapsis
const satEllipses = { A: toUnits(EARTH_RADIUS_KM + 417), EC: 0.0007 };
const satOrbitPath = getOrbitPath(satEllipses.A, satEllipses.EC);
satelliteAbsidalOrbitalPlane.add(satOrbitPath);


// Add additional elements if debug is enabled
if (DEBUG)
    init_debug(scene, tiltedEarth, earth, satellite);

// Initial/default camera view
const defaultCameraPos = (function() {
  const cameraDistance = toUnits(1.2 * MOON_DISTANCE_KM); // From Earth center
  const cameraAzimuth = toRadians(90);  // Degrees clockwise from +Z axis
  const cameraElevation = toRadians(15); // Degrees above horizontal plane
  return new THREE.Vector3(
    cameraDistance * Math.cos(cameraElevation) * Math.sin(cameraAzimuth),
    cameraDistance * Math.sin(cameraElevation),
    cameraDistance * Math.cos(cameraElevation) * Math.cos(cameraAzimuth)
  );
})();
// Here we only instantiate the view mananger, we must init the main view later
// when we know Earth position
const views = new ViewManager(scene, renderer, defaultCameraPos, earth);

// Init observer object
const observer = new Observer(views, false);

// Init satObserver object
const satObserver = new Observer(views, true);

// List of objects that can have a locked view on them
const trackableBodies = [earth, moon];

// Return any trackable object under mouse
function pickObject(mouseX, mouseY) {

  const raycaster = new THREE.Raycaster();
  const mouse = new THREE.Vector2();

  // Convert mouse position to normalized device coordinates (-1 to +1)
  mouse.x =  (mouseX / window.innerWidth)  * 2 - 1;
  mouse.y = -(mouseY / window.innerHeight) * 2 + 1;

  // Update the raycaster with the camera and mouse position
  const view = views.getActive();
  raycaster.setFromCamera(mouse, view.camera);

  // Find intersections with candidates
  const intersects = raycaster.intersectObjects(trackableBodies);
  if (intersects.length === 0)
    return null;

  return {
    object: intersects[0].object,
    point: intersects[0].point // point on the surface in world space
  };
}

// Lock the camera into a geostationary orbit around the selected object
function lockToOrbit(objectToLock, surfaceWorldPoint) {

  // Set a marker in local coordinates on object's surface point
  const marker = getMarker(toUnits(50), 0x00ff00, surfaceWorldPoint);
  objectToLock.worldToLocal(marker.position);
  objectToLock.add(marker);

  // Lock the active view aligned to object's center
  views.lockToOrbitView(marker, objectToLock);
}

// ============================================================================
// CELESTIAL MECHANICS
// ============================================================================

// Convert a position vector from the Scene's World Frame back to the
// standard ICRS Ecliptic J2000 Frame and rescale from units to km.
function convertToJPL(worldPosition) {
  const x =   worldPosition.x; // x_icrs =  x_scene
  const y = - worldPosition.z; // y_icrs = -z_scene
  const z =   worldPosition.y; // z_icrs =  y_scene
  return [fromUnits(x), fromUnits(y), fromUnits(z)];
}

// Convert an offset vector in world frame into a position in a
// local frame with a given origin
function offsetToLocal(offset, origin) {
  origin.getWorldPosition(temp1Vec);
  temp1Vec.add(offset);
  origin.worldToLocal(temp1Vec);
  offset.copy(temp1Vec);
}

// Distance from focus given true anomaly (ellipse polar equation)
function radius(ta, A, EC) {
  return A * (1 - EC * EC) / (1 + EC * Math.cos(ta));
}

// Computes an array of 3D vertex on an elliptic curve in the Moon's absidal
// XZ orbital plane, with origin at focus and +X axis pointing toward periapsis
function ellipticCurve(A, EC) {
  const SEGMENTS = 360;
  const points = new Float32Array(SEGMENTS * 3);

  // Loop starts at periapsis and moves in counter-clockwise direction
  // at i = 0 -> [r, 0, 0]
  for (let i = 0; i < SEGMENTS; i++) {
    const ta = (i / SEGMENTS) * 2 * Math.PI; // true anomaly
    const r = radius(ta, A, EC); // distance from focus at true anomaly
    points[i*3 + 0] = r * Math.cos(ta); // x
    points[i*3 + 1] = 0;                // y  (on orbital plane)
    points[i*3 + 2] =-r * Math.sin(ta); //-z due to local frame conventions
  }
  return points;
}

// Solves Kepler's Equation for Eccentric Anomaly (E) with
// the Newton-Raphson iterative method
function solveKepler(M, e) {

  // Kepler's Equation is: M = E - e * sin(E)
  const maxIterations = 10;
  const tolerance = 1e-6; // A small tolerance for accuracy
  let E = M; // Initial guess for E is M

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

// Solve for position in a 2D XZ orbital plane, with origin
// at focus and +X axis pointing toward periapsis
function solveAbsidalOrbit(current) {

  // Solve Kepler's equation for Eccentric Anomaly (E)
  const M = current.M; // Mean anomaly
  const EC = current.EC; // Eccentricity
  const E = solveKepler(M, EC);

  // Calculate True Anomaly (TA), the actual angle from periapsis
  const [sinE, cosE] = [Math.sin(E), Math.cos(E)];
  const TA = Math.atan2(
    Math.sqrt(1 - EC * EC) * sinE, // z-coordinate on auxiliary circle
    cosE - EC                      // x-coordinate on auxiliary circle
  );

  // Radius as distance from focus
  const r = current.A * (1 - current.EC * cosE);

  // Calculate the 2D position in the XZ orbital plane
  // The plane's +X axis is oriented toward periapsis so we use TA
  const x =  r * Math.cos(TA);
  const z = -r * Math.sin(TA); // -z for prograde motion
  return [x, z];
}

// Calculate the EMB/Earth position in its absidal orbital plane,
// with +X axis pointing toward periapsis / perihelion
function getEBMSunPosition(current) {

  // Calculate Mean Anomaly (M) from Mean Longitude from Vernal Equinox (L)
  // and Longitude of Perihelion from VE (W_bar)
  current.M = current.L - current.W_bar;

  const [x, z] = solveAbsidalOrbit(current);

  // Return 3D coordinates in the local plane's frame
  return [x, 0, z];
}

// Calculate the Moon position in the Moon absidal orbital plane,
// with +X axis pointing toward periapsis / perigee
function getMoonEBMPosition(current) {

  // Mean Anomaly is already updated
  current.M = current.MA;

  const [x, z] = solveAbsidalOrbit(current);

  // Return 3D coordinates in the local plane's frame
  return [x, 0, z];
}

// Calculate the satellite position in its absidal orbital plane,
// with +X axis pointing toward periapsis / perigee
function getSatellitePosition(current) {

  // Mean Anomaly is already updated
  current.M = current.MA;

  const [x, z] = solveAbsidalOrbit(current);

  // Return 3D coordinates in the local plane's frame
  return [x, 0, z];
}

// Return the height above visible horizon of a give target
function heightAboveHorizon(target, radius) {

  // Get target's position in the observer frame
  target.getWorldPosition(temp1Vec);
  const targetLocalPos = observer.marker.worldToLocal(temp1Vec);

  // Correction for the object's upper limb (its radius)
  const radiusCorrection = radius;

  // Correction for atmospheric refraction (lifts the image)
  let refractionCorrection = 0;
  if (observer.object === earth) {
    const distance = targetLocalPos.length();
    refractionCorrection = distance * Math.tan(HORIZON_REFRACTION);
  }

  // The object is visible if its center is above this negative threshold
  const visibilityThreshold = -(radiusCorrection + refractionCorrection);
  return targetLocalPos.y - visibilityThreshold;
}


// ============================================================================
// ANIMATION LOOP
// ============================================================================

// Clock works also when tab is hidden, can be set/reset by UI, time
// can run faster and/or backward from real time and is adjustable
const SimClock = class {
  constructor(masterEpoch) {
    this.masterEpochTime = masterEpoch.getTime(); // our t=0
    this.reset();
    this.running = true;
  }
  reset() {
    this.simulationTimeAtPause = Date.now();
    this.realTimeAtResume = this.simulationTimeAtPause;
    this.speedX = SIM_TIME_SPEED_UP;
  }
  setDate(newDate) {
    // Adjust the "banked" time to make the total time equal the newDate
    this.simulationTimeAtPause = newDate.getTime() - this._deltaSinceLastResume();
  }
  _deltaSinceLastResume() {
    // Accumulate time only while running. It is simulation time, not real time
    return this.running ? this.speedX * (Date.now() - this.realTimeAtResume) : 0;
  }
  getTime() {
    const simTime = this.simulationTimeAtPause + this._deltaSinceLastResume();
    return simTime; // in msecs
  }
  elapsed() {
    // Initial offsets in simulation are computed on master epoch
    return this.getTime() - this.masterEpochTime;
  }
  togglePause() {
    if (this.running) {
      // Before stopping, add the "ticking" time to the "banked" time
      this.simulationTimeAtPause += this._deltaSinceLastResume();
    } else {
      this.realTimeAtResume = Date.now();
    }
    this.running = !this.running;
    return this.running;
  }
  speed() {
    return this.speedX / SIM_TIME_SPEED_UP; // SIM_TIME_SPEED_UP -> speed == 1
  }
  setSpeed(speed) {
    // Bank any accumulated time at the old speed by simulating a
    // stop/start, then set the new speed. It works whether the
    // clock is running or paused.
    this.togglePause();
    this.togglePause();
    this.speedX = SIM_TIME_SPEED_UP * speed;
  }
};

// Main animation loop
function animate(simulation) {

  // Our simulation output, apart from rendering
  let simStepData = null;

  // Three.js rotation order is 'XYZ' but astronomical standard defines
  // the orientation of an orbital plane in a specific sequence:
  // 1. Rotate around the Pole
  // 2. Tilt around the Line of Nodes

  // Elapsed seconds and julian centuries since J2000 Epoch
  const elpasedMsec = simulation.clock.elapsed();
  const elapsedSeconds = elpasedMsec / 1000;
  const julianCenturies = elapsedSeconds / (100 * JULIAN_YEAR);

  // Propagate orbital parameters to current simulation time
  const currentEarth = propagateEarth(julianCenturies);
  const currentMoon = propagateMoon(julianCenturies);

  // Initial static rotations of the Sun orbital plane
  // Set the Euler rotation order to 'YXZ' to match the
  // astronomical convention for an orbital plane.
  earthOrbitalPlane.rotation.order = 'YXZ';
  earthOrbitalPlane.rotation.x = currentEarth.IN;
  earthOrbitalPlane.rotation.y = currentEarth.W_bar;

  // Compute EMB position in the Earth-Sun orbital plane
  const [xe, ye, ze] = getEBMSunPosition(currentEarth);
  embPosWorldVec.set(xe, ye, ze);

  // Use a rescaled EMB-Sun distance for rendering purposes
  if (!simulation.validateMode) {
    embPosWorldVec.normalize().multiplyScalar(SUN_LIGHT_DISTANCE);
  }
  embPivot.position.copy(embPosWorldVec);

  // Earth (apsidal) orbital plane has +X axis pointing toward
  // perihelion but we want EMB to be in world frame, because
  // Earth's axial tilt orientation depends on it, so invert
  // the parent's quaternion
  earthOrbitalPlane.updateWorldMatrix(true, false);
  earthOrbitalPlane.getWorldQuaternion(toWorldQuat);
  embPivot.quaternion.copy(toWorldQuat.invert());

  // Inclination to Ecliptic and ascending node precession
  // in backwards (retrograde) direction. Reference frame
  // is Heliocentric Ecliptic at J2000 Epoch
  moonNodesOrbitalPlane.rotation.order = 'YXZ';
  moonNodesOrbitalPlane.rotation.x = currentMoon.IN;
  moonNodesOrbitalPlane.rotation.y = currentMoon.OM;

  // Set origin of Moon's nodes orbital plane onto EMB, this is
  // only a translation, local frame axis directions do not change
  embPivot.getWorldPosition(embPosWorldVec);
  moonNodesOrbitalPlane.position.copy(embPosWorldVec);

  // Ascending and descending nodes lie on the x-axis of the parent
  // orbital plane of nodes and due to Moon's absidal precession
  // their distance to Earth changes dynamically.
  const ta = -currentMoon.W; // - argument of perigee by definition
  const rMoonAsc = radius(ta, currentMoon.A, currentMoon.EC);
  const rMoonDes = radius(ta + Math.PI, currentMoon.A, currentMoon.EC);
  ascendingNodeMarker.position.set(rMoonAsc, 0, 0);
  descendingNodeMarker.position.set(-rMoonDes, 0, 0);

  // Absidal precession in backwards (retrograde) direction
  moonAbsidalOrbitalPlane.rotation.y = currentMoon.W;

  // Compute Moon position in the Moon's absidal plane
  const [xm, ym, zm] = getMoonEBMPosition(currentMoon);
  tiltedMoon.position.set(xm, ym, zm);

  // Moon's axis tilt
  // Cassini's 3rd Law  dictates that the Moon's spin axis, its orbital pole,
  // and the Ecliptic Pole are always coplanar, so pole of ecliptic (~1.54°)
  // is always *between* the 2 poles of orbit (~5.1°) and spin (~6.68°) and
  // the Moon's spin axis is "tilted back" from its orbital plane. Moon orbital
  // plane inclination to Ecliptic is done rotating around its local X axes of
  // a positive angle, so we rotate around X axis with an opposite sign.
  tiltedMoon.rotation.x = - toRadians(MOON_AXIAL_TILT);

  // Moon rotation around its Pole axis, in tidal locking with Earth
  // Moon texture map is loaded centered on 0° longitude, the Prime Meridian of
  // the Moon is aligned with its local +X axis.
  // Moon's absidal local frame +X axis direction goes from the EMB to the Perigee,
  // so when the Moon is at the Perigee the 0° longitude lies on the far side
  // of the Moon. We rotate around Y axis of 180 degrees to move it facing Earth.
  const moonRotationSpeed = 2 * Math.PI / MOON_ORBITAL_PERIOD;
  moon.rotation.y = Math.PI + moonRotationSpeed * elapsedSeconds;

  // Earth's axial tilt and retrograde (clockwise) precession.
  // Tilt is a negative rotation on the X-axis (Vernal Point direction). This orients
  // the North Pole towards the Winter Solstice direction, correctly setting up
  // the seasons. Precession is a slow, negative rotation around the tilted Y-axis.
  const earthPrecessionSpeed = 2 * Math.PI / EARTH_AXIAL_PERIOD;
  tiltedEarth.rotation.x = - toRadians(EARTH_AXIAL_TILT);
  tiltedEarth.rotation.y = - earthPrecessionSpeed * elapsedSeconds;

  // Earth daily rotation
  // Earth texture map is loaded centered on 0° longitude, the Prime Meridian of
  // the Earth is aligned with its local +X axis that is also the world +X axis
  // and points from Earth toward Vernal Point. Greenwich Apparent Sidereal
  // Time (GAST) is the angle between the Vernal Point and the Prime Meridian
  // at the GAST epoch time.
  const earthRotationSpeed = 2 * Math.PI / SIDERAL_DAY;
  const timseSinceGAST = elapsedSeconds - GAST_TIME; // in secs
  earth.rotation.y = GAST + earthRotationSpeed * timseSinceGAST;

  // Earth placement according to barycenter rule
  //
  // Earth lies on the Moon's orbital plane, where EMB is at the Origin,
  // so we rescale and invert the EBM->Moon vector and transform in
  // local coord, this is the Earth's new position in its local frame
  tiltedMoon.getWorldPosition(moonPosWorldVec);
  embPivot.getWorldPosition(embPosWorldVec);
  earthPosWorldVec.subVectors(moonPosWorldVec, embPosWorldVec) // EMB->Moon
                  .multiplyScalar(MOON_EARTH_MASS_RATIO)
                  .negate(); // EMB->Earth (in world frame)
  offsetToLocal(earthPosWorldVec, embPivot);
  tiltedEarth.position.copy(earthPosWorldVec);

  // Satellite placement
  //
  const timseSinceSatEpoch = elapsedSeconds - ISS_TIME; // in secs
  const currentSatellite = propagateSatellite(ISS_J2000, timseSinceSatEpoch);

  // Apply orientation to the satellite's orbital planes
  satelliteNodesOrbitalPlane.rotation.order = 'YXZ';
  satelliteNodesOrbitalPlane.rotation.x = currentSatellite.IN;
  satelliteNodesOrbitalPlane.rotation.y = currentSatellite.OM;
  satelliteAbsidalOrbitalPlane.rotation.y = currentSatellite.W;

  // Compute Satellite position in its absidal plane
  const [xs, ys, zs] = getSatellitePosition(currentSatellite);
  satellite.position.set(xs, ys, zs);

  // Force hierarchy recalculation and recompute real positions before
  // second part of animation
  scene.updateWorldMatrix(true, true);

  if (simulation.validateMode) {
    const m = tiltedMoon.getWorldPosition(moonPosWorldVec);
    const e = tiltedEarth.getWorldPosition(earthPosWorldVec);
    const emb = embPivot.getWorldPosition(embPosWorldVec);

    // Translate frame origin to EMB
    const em = m.clone().sub(emb); // EMB -> Moon (in world frame)
    //offsetToLocal(em, embPivot); // EMB -> Moon (in EMB frame)

    const moonV  = convertToJPL(em);
    const earthV = convertToJPL(e);
    const timeV = elpasedMsec + simulation.clock.masterEpochTime; // UTC time
    return {
      time: timeV,         // msec
      sunPosition: earthV, // km
      moonPosition: moonV  // km
    };
  }

  // Set Sun light to look at Earth
  tiltedEarth.getWorldPosition(earthPosWorldVec);
  sunLight.target.position.copy(earthPosWorldVec);

  // Update the end point of the line to match the sun's new position
  if (DEBUG) {
    set_sunline_length(earthPosWorldVec);
  }

  // In dry run mode sample the proper value and return
  if (simulation.dryRunFunction) {
      return simulation.dryRunFunction();
  }

  // Update celestial events
  const isObserver = (observer.object !== null);
  const timeForward = (simulation.clock.speed() > 0);
  simulation.eventManager.update(timeForward, isObserver);

  // If we are in observer or satellite view get
  // specific data
  if (simulation.isObserverView()) {
    simStepData = observer.getGeoData()
  } else if (simulation.isSatelliteView()) {
    simStepData = satObserver.getSatData(elapsedSeconds);
  }

  // One time event to reset camera view after a change of date
  if (simulation.view_needs_reset) {
    const view = views.getActive();
    view.reset(earthPosWorldVec);
    simulation.view_needs_reset = false;
  }

  // If is first frame init main view and other stuff that requires
  // updated orbital parameters and Earth position
  if (views.getActive() === null) {
    views.init(earthPosWorldVec);
  }

  // Update the active controls and render with the active camera
  const activeCamera = views.update();
  renderer.render(scene, activeCamera);

  return simStepData;
}

// ============================================================================
// EVENT DETECTORS
// ============================================================================

// Functions used in dry run to sample the requested simulation values
function sunHeight() {
  // Calculate the equivalent radius using the captured real distance
  const sunDistance = earthPosWorldVec.length();
  const equivalentSunRadius = toUnits(SUN_RADIUS_KM) * (SUN_LIGHT_DISTANCE / sunDistance);
  return heightAboveHorizon(sunLight, equivalentSunRadius);
}

function moonHeight() {
  return heightAboveHorizon(tiltedMoon, toUnits(MOON_RADIUS_KM));
}

// Calculates the signed distance of a target object from a reference plane
function signedDistanceToPlane(target, planeAxisObject) {

  // Use the local Y-axis of the object as the defining axis
  temp1Vec.set(0, 1, 0);
  const planeYAxis = temp1Vec.transformDirection(planeAxisObject.matrixWorld);
  const targetPosWorldVec = target.getWorldPosition(temp2Vec);

  // The reference plane's normal is perpendicular to both the
  // Sun-Earth vector and the defining axis
  const planeNormal = temp3Vec
      .crossVectors(earthPosWorldVec, planeYAxis)
      .normalize();

  // The dot product gives the perpendicular distance from the target to the plane
  const signedDistanceFormPlane = planeNormal.dot(targetPosWorldVec);
  return signedDistanceFormPlane;
}

// Calculates the distance from the observer's local meridian.
// A zero-crossing indicates a solar transit (noon or midnight).
// The plane is defined by Earth's rotational axis
function sunTransitDistance() {
  return signedDistanceToPlane(observer.marker, earth);
}

// Calculates the distance of the Moon from the Sun-Earth line, as projected
// onto the ecliptic plane. A zero-crossing indicates a syzygy (New or Full Moon).
// The plane is defined by the ecliptic pole (the normal to Earth's orbital plane).
function newMoonDistance() {
  return signedDistanceToPlane(moon, earthOrbitalPlane);
}

// ============================================================================
// EXPORT CLASS
// ============================================================================

class Simulation {
  constructor(validateMode) {
    // Don't render and return Moon / Earth positions
    this.validateMode = validateMode;
    this.dryRunFunction = null;
    const eventFuns = [sunHeight, moonHeight, sunTransitDistance, newMoonDistance];

    // Our simulation clock, init with J2000 Epoch
    this.clock = new SimClock(J2000_EPOCH);

    // Create a new event manager and expose its 'on' method
    this.eventManager = new CelestialEventManager(this, eventFuns);
    this.on = this.eventManager.on.bind(this.eventManager);

    // A flag to handle deferred view update after date change
    this.view_needs_reset = false;

    this.update = () => animate(this);
    this.getTime = this.clock.getTime.bind(this.clock);
    this.speed = this.clock.speed.bind(this.clock);
    this.setSpeed = this.clock.setSpeed.bind(this.clock);
    this.togglePause = this.clock.togglePause.bind(this.clock);
    this.disposeView = views.dispose.bind(views);
    this.lockToOrbit = lockToOrbit;
    this.pickObject = pickObject;
  }
  setDryRunFunction(fun) {
    this.dryRunFunction = fun;
  }
  getRenderer() {
    return renderer;
  }
  reverseSpeed() {
    const speed = -this.clock.speed();
    this.clock.setSpeed(speed);
    return speed;
  }
  reset() {
    this.clock.reset();
    this.setActiveView(0);
  }
  saveView() {
    // Save camera position to keep current view
    this.eventManager.reset();
    tiltedEarth.getWorldPosition(earthPosWorldVec);
    const view = views.getActive();
    view.setDefault(earthPosWorldVec);
  }
  setDate(newDate) {
    // While in dry run we don't want to override any pending reset
    if (!this.view_needs_reset) {
      this.saveView();
      this.view_needs_reset = true; // deferred to after simulation step
    }
    this.clock.setDate(newDate);
  }
  goToNextEvent(eventName, goNext, timeForward) {
    this.saveView();
    this.view_needs_reset = true; // deferred to after simulation step
    this.eventManager.goToNextEvent(eventName, goNext, timeForward);
  }
  setActiveView(viewIndex) {
    views.setActive(viewIndex);
    // Hide orbit path clutter and set realistic scale for satellite
    const isLocalView = this.isSatelliteView() || this.isObserverView();
    satOrbitPath.visible = !isLocalView;
    const scale = isLocalView ? SAT_REALISTIC_SCALE : 1.0;
    satellite.scale.set(scale, scale, scale);
  }
  enterSatelliteView() {
    const eyeHeight = toUnits(0);   // km above surface
    const lookAhead = toUnits(100); // look at km ahead on horizon

    const viewIndex = satObserver.createObserverView(earth, satellite, eyeHeight, lookAhead);
    this.setActiveView(viewIndex);
  }
  exitSatelliteView() {
    // FIXME actively change view to last one before dispose
    satObserver.disposeObserverView();
  }
  isSatelliteView() {
    const view = views.getActive();
    return view !== null && view === satObserver.getView();
  }
  createObserverView(object, surfaceWorldPoint) {
    const eyeHeight = toUnits(5); // km above surface
    const lookAhead = toUnits(100); // look at km ahead on horizon
    this.eventManager.reset();
    const marker = getMarker(toUnits(10), 0xFFFFFF, surfaceWorldPoint);
    object.worldToLocal(marker.position);
    object.add(marker);
    const viewIndex = observer.createObserverView(object, marker, eyeHeight, lookAhead);
    return viewIndex;
  }
  placeObserverAt(latDeg, lonDeg) {
    this.eventManager.reset();
    observer.placeAt(latDeg, lonDeg);
  }
  disposeObserverView() {
    this.eventManager.reset();
    const marker = observer.marker;
    marker.removeFromParent();
    disposeMarker(marker);
    observer.disposeObserverView();
  }
  isObserverView() {
    const view = views.getActive();
    return view !== null && view === observer.getView();
  }
  isOrbitLocked(object) {
    const lockedObjects = views.getOrbitLockedObjects();
    return lockedObjects.includes(object);
  }
  unlockCamera() {
    const view = views.getActive();
    const marker = view.getLockedObject();
    marker.removeFromParent();
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
