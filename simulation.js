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
import { init_debug, set_sunline_length, scene_axes } from './validate.js';

export const DEBUG = false; // Enable axes and objects visualizations for debug purposes
export const VALIDATE = false; // Run validation script instead of normal visualization mode

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
const SUN_RADIUS_KM = 695700;
const MOON_DISTANCE_KM = 384400;    // average
const EARTH_AXIAL_TILT = 23.439291; // degrees
const MOON_AXIAL_TILT = 6.68;       // degrees

// Earth-Moon Barycenter (EBM) Correction Constants to convert
// the EMB->Moon distance to the EMB->Earth distance
//  m₁ × d₁ = m₂ × d₂  --> EBM_E = EBM_M * MOON_EARTH_MASS_RATIO
const MOON_EARTH_MASS_RATIO = 1 / 81.3;

// The official IAU 2012 definition of the Astronomical Unit in kilometers.
const AU = 149597870.7;

// Scale conversions
const KM_PER_UNIT = EARTH_RADIUS_KM / 50; // Earth radius to units
const toUnits = km => km / KM_PER_UNIT;
const fromUnits = u => u * KM_PER_UNIT;
const toRadians = degrees => degrees * Math.PI / 180;
const toDegrees = rad => rad * 180 / Math.PI;
const HMSToRadians = s => s.split(/\s+/).reduce((acc, v, i) => acc + parseFloat(v) / [1, 60, 3600][i], 0) * (Math.PI / 12);

// Rescaled Sun distance for placing directional Sun light source
const SUN_LIGHT_DISTANCE = toUnits(10 * MOON_DISTANCE_KM);

// Time constants for simulation
const SIM_TIME_SPEED_UP = 60; // 1 simulation hour elapses in 1 minute

// Master Epoch of JPL data sources
const MASTER_EPOCH = new Date('2025-09-10T09:34:03Z');
const J2000_EPOCH = new Date('2000-01-01T12:00:00Z');
const VERNAL_EQUINOX_2025 = new Date('2025-03-20T09:01:00Z');
const VERNAL_EQUINOX_JPL = new Date('2000-03-20T07:26:28Z'); // 2000-Mar-20 07:26:28 TDB

const HOUR = 3600; // In SI seconds
const SOLAR_DAY = 24 * HOUR; // Mean Solar Day
const JULIAN_YEAR = 365.25 * SOLAR_DAY;

const SIDERAL_DAY = 23.9344696 * HOUR; // Period of a 360 degrees rotation
const SIDEREAL_YEAR = 365.256363004 * SOLAR_DAY; // For epoch J2000.0

const EARTH_AXIAL_PERIOD = 25772 * JULIAN_YEAR;

const MOON_ORBITAL_PERIOD = 27.32166 * SOLAR_DAY;
const MOON_NODE_PERIOD = 18.612815932 * JULIAN_YEAR;
const MOON_ABSIDAL_PERIOD = 8.85 * JULIAN_YEAR;

// Calculate speed constants in radians/sec
const EARTH_ROTATION_SPEED = 2 * Math.PI / SIDERAL_DAY;
const MOON_ROTATION_SPEED = 2 * Math.PI / MOON_ORBITAL_PERIOD; // Tidal lock
const EARTH_AXIAL_SPEED = 2 * Math.PI / EARTH_AXIAL_PERIOD;

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

// Target Body: Moon, Coordinate Center: Earth-Moon Barycenter [EMB]
const MOON_J2000 = {
  A: toUnits(379700),    // Semi-major axis in Km
 EC: 0.0554,             // Eccentricity
 IN: toRadians(5.145),   // Inclination to Ecliptic
 OM: toRadians(125.044), // Longitude of the Ascending Node (from VE)
  W: toRadians(318.15),  // Argument of periapsis (from AN)
 MA: toRadians(134.963), // Mean anomaly
};

// Rate of Change (per Julian century)
const MOON_CHANGE_RATES = {
  A: 0,
 EC: 0,
 IN: 0,
 OM:-toRadians(100 * JULIAN_YEAR * 360 / MOON_NODE_PERIOD),    // clockwies (moves backwards)
  W: toRadians(100 * JULIAN_YEAR * 360 / MOON_ABSIDAL_PERIOD), // counter-clockwise (moves forwards)
 MA: toRadians(100 * JULIAN_YEAR * 360 / MOON_ORBITAL_PERIOD), // counter-clockwise

 // const M = 134.96292 + 477198.86753 * T + 0.00924 * T * T;
}

/*
    Moon Ephemeris (osculating elements) from JPL Horizons

    https://ssd.jpl.nasa.gov/horizons/app.html

    The reference frame of data source is the ICRS Ecliptic at J2000 Epoch
    with origin at Earth-Moon Barycenter (EBM)

    EPOCH: JPL osculating elements refer to MASTER_EPOCH time

    To solve a 2D orbital ellipses we use a 3D local frame defined as following:

    ORIENTATION: Local frame is a right-handed coordinate system

    ORIGIN: The origin (0,0,0) is the nearest focus of the orbital ellipses

    THE ORBITAL PLANE: The local XZ plane

    AXES DEFINITION:
      +X Axis: This is the principal direction, defined as the "line of nodes".
               It points from the focus towards the Ascending Node

      +Y Axis: Perpendicular to the orbital plane, this is the UP direction
               pointing toward orbital north pole

      +Z Axis: At the Ascending Node a counter-clockwise orbit (seen from +Y)
               has velocity with direction -Z
*/

// Target Body: Moon [Luna], Coordinate Center: Earth-Moon Barycenter [500@3]
const MOON_DATA = {
     A: toUnits(3.796241532734630E+05),
    EC: 5.078618152688460E-02,
    IN: toRadians(5.280291370607761E+00),
    OM: toRadians(3.479716750341908E+02),
     W: toRadians(3.711163694612181E+01),
    MA: toRadians(1.290172290680006E+01),
    PR: 2.343381844716793E+06,
};

// Target Body: Earth, Coordinate Center: Earth-Moon Barycenter [500@3]
const EARTH_DATA = {
    //  A: toUnits(4.669391144219639E+03),
    // EC: 5.078618152693617E-02,
    // IN: toRadians(5.280291370607768E+00),
    // OM: toRadians(3.479716750341909E+02),
    //  W: toRadians(2.171116369461202E+02),
    // MA: toRadians(1.290172290678550E+01),
    // PR: 2.343381844717038E+06,

    L_Ap_Sid_Time: HMSToRadians("08 52 45.4526"), // Greenwich Apparent Sidereal Time (GAST)
};

// ============================================================================
// SCENE SETUP
// ============================================================================

/*
  Setup a hierarchy of 3D Objects to split complex movement in independent
  rotations
*/

// Helper to create a marker object
function getMarker(color, position) {
  const geometry = new THREE.SphereGeometry(toUnits(600));
  const material = new THREE.MeshBasicMaterial({ color: color });
  const marker = new THREE.Mesh(geometry, material);
  marker.position.set(...position);
  return marker;
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

// Delayed init after inital orbital parameters have been updated
function init_moon_path(current) {
  // An elliptical curve to show Moon's orbit, local +X axis points
  // toward periapsis / perigee
  const [points, rMoonAsc, rMoonDes] = ellipticCurve(current);
  const vertexBuffer = new THREE.BufferAttribute(points, 3);
  const orbitPathGeo = new THREE.BufferGeometry().setAttribute('position', vertexBuffer);
  const orbitPathMat = new THREE.LineBasicMaterial({
    color: 0x4488ff, // light blue
    transparent: true,
    opacity: 0.35
  });
  const moonOrbitPath = new THREE.LineLoop(orbitPathGeo, orbitPathMat);
  moonAbsidalOrbitalPlane.add(moonOrbitPath);

  // Add markers on the Moon's orbit ascending and descending nodes
  moonNodesOrbitalPlane.add(getMarker(0xff00ff, [+rMoonAsc, 0, 0])); // Magenta
  moonNodesOrbitalPlane.add(getMarker(0xffff00, [-rMoonDes, 0, 0])); // Yellow
}

// Add additional elements if debug is enabled
if (DEBUG)
    init_debug(scene, tiltedEarth, earth);

// Initial/default camera view
const defaultCameraPos = (function() {
  const cameraDistance = toUnits(2 * MOON_DISTANCE_KM); // From Earth center
  const cameraAzimuth = toRadians(180);  // Degrees clockwise from +Z axis
  const cameraElevation = toRadians(10); // Degrees above horizontal plane
  return new THREE.Vector3(
    cameraDistance * Math.cos(cameraElevation) * Math.sin(cameraAzimuth),
    cameraDistance * Math.sin(cameraElevation),
    cameraDistance * Math.cos(cameraElevation) * Math.cos(cameraAzimuth)
  );
})();
// Here we only instantiate the view maanger, we must init the main view later
// when we know Earth position
const views = new ViewManager(scene, renderer, defaultCameraPos);

// ============================================================================
// OBSERVER VIEW MANAGEMENT
// ============================================================================

// List of objects that can have a locked view on them
const trackableBodies = [earth, moon];

// Tracks state for an observer on a planet
const observerState = {
    object: null,
    marker:null, // in local coordinates
    observerOnEarth: false,
    viewIndex: null,
    tempVec: new THREE.Vector3(),
};

// Helper used by observer and orbit locking
function createMarker(radius, color) {
  return new THREE.Mesh(
    new THREE.SphereGeometry(radius),
    new THREE.MeshBasicMaterial({ color })
  );
}

// Release marker resources
function disposeMarker(marker) {
  marker.geometry.dispose();
  marker.material.dispose();
}

// Add a child object, usually a marker, on a surface of
// a body and align local coordinates of added child
function attachToSurface(marker, object, surfacePoint) {

  // Surface point must be in world coords
  marker.position.copy(surfacePoint);
  object.worldToLocal(marker.position);
  orientToSurface(marker);
  if (marker.parent !== object) // in case we just move the marker
    object.add(marker);
  marker.updateWorldMatrix(true, false);
}

// Orient the local coordinate frame of an object on the surface of a parent body
function orientToSurface(marker) {
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
function lockToOrbit(objectToLock, surfacePoint) {

  // Set a marker in local coordinates on object's surface point
  const marker = createMarker(toUnits(50), 0x00ff00);
  attachToSurface(marker, objectToLock, surfacePoint);

  // Lock the active view aligned to object's center
  views.lockToOrbitView(marker, objectToLock);
}

// Place an observer on an object surface
function enterObserverMode(object, surfacePoint) {

  // Currently we handle only one observer
  console.assert(observerState.object === null, "Setting already existing observer");

  // Position and orient the marker on the surface
  const marker = createMarker(toUnits(10), 0xFFFFFF);
  attachToSurface(marker, object, surfacePoint);

  observerState.object = object;
  observerState.marker = marker;

  // Init observer events
  observerState.observerOnEarth = (object === earth);

  // Create a new observer view placed on marker
  const eyeHeight = toUnits(5); // km above surface
  const lookAhead = toUnits(100); // look at km ahead on horizon
  const viewIndex = views.createObserverView(marker, eyeHeight, lookAhead);
  observerState.viewIndex = viewIndex;

  // Lock the view camera on the marker
  const view = views.get(viewIndex);
  view.lockTo(marker, true);

  return viewIndex;
}

// Move the marker to a new position on planet surface
function placeObserverAt(latDeg, lonDeg) {

  const view = views.get(observerState.viewIndex);
  if (views.getActive() !== view)
    return;

  const latRad = toRadians(latDeg);
  const lonRad = toRadians(lonDeg);
  const radius = observerState.marker.position.length();

  // Convert spherical coordinates to Cartesian (x, y, z) in the planet's local
  // frame. Y is the polar axis. +X is the prime meridian (0° longitude).
  // The -z for sine is to match the longitude calculation in getGeoData.
  const newSurfacePoint = new THREE.Vector3();
  newSurfacePoint.y = radius * Math.sin(latRad);
  const xzRadius = radius * Math.cos(latRad); // Radius of the circle at this latitude
  newSurfacePoint.x = xzRadius * Math.cos(lonRad);
  newSurfacePoint.z = -xzRadius * Math.sin(lonRad);

  // Move the marker to new surface point and realign its plane
  // Surface point must be in world coords
  const object = observerState.object;
  object.localToWorld(newSurfacePoint)
  attachToSurface(observerState.marker, object, newSurfacePoint);

  // Camera is already locked to the marker, locking update will
  // set correct camera and target position at next frame.
}

// Drop observer view and release relative resources
function exitObserverMode() {
  const view = views.get(observerState.viewIndex);
  view.unlock();
  const marker = observerState.marker;
  marker.parent.remove(marker);
  disposeMarker(marker);
  views.dispose(observerState.viewIndex);
  observerState.object = null;
  observerState.viewIndex = null;
  observerState.marker = null;
}

// ============================================================================
// CELESTIAL MECHANICS
// ============================================================================

const tempVec = new THREE.Vector3();

// Convert a position vector from world frame to the JPL Ecliptic J2000 frame
// and from units values back to km
function convertToJPL(position) {

  const jpl_x =   position.x; // Matches directly
  const jpl_y = - position.z; // Negate the Z component
  const jpl_z =   position.y; // Y is the UP axis

  return [fromUnits(jpl_x), fromUnits(jpl_y), fromUnits(jpl_z)];
}

// Convert an offset vector in world frame into a position in a
// local frame with a given origin
function offsetToLocal(offset, origin) {
  origin.getWorldPosition(tempVec);
  tempVec.add(offset);
  origin.worldToLocal(tempVec);
  offset.copy(tempVec);
}

// Computes an array of 3D vertex on an elliptic curve in the Moon's absidal
// XZ orbital plane, with origin at focus and +X axis pointing toward periapsis
function ellipticCurve(current) {
  const SEGMENTS = 360;
  const points = new Float32Array(SEGMENTS * 3);
  const A = current.A;
  const EC = current.EC;
  const W = current.W;

  // Distance from focus given true anomaly (ellipse polar equation)
  const radius = ta => A * (1 - EC * EC) / (1 + EC * Math.cos(ta));

  // Loop starts at periapsis and moves in counter-clockwise direction
  // at i = 0 -> [r, 0, 0]
  for (let i = 0; i < SEGMENTS; i++) {
    const ta = (i / SEGMENTS) * 2 * Math.PI; // true anomaly
    const r = radius(ta); // distance from focus at true anomaly
    points[i*3 + 0] = r * Math.cos(ta); // x
    points[i*3 + 1] = 0;                // y  (on orbital plane)
    points[i*3 + 2] =-r * Math.sin(ta); //-z due to local frame conventions
  }

  // Distances form origin of ascending and descending nodes
  const rAsc = radius(-W);
  const rDes = radius(-W + Math.PI);
  return [points, rAsc, rDes];
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

function heightAboveHorizon(target, radius) {

  // Get target's position in the observer frame
  const temp = observerState.tempVec;
  target.getWorldPosition(temp);
  const targetLocalPos = observerState.marker.worldToLocal(temp);

  // Correction for the object's upper limb (its radius)
  const radiusCorrection = radius;

  // Correction for atmospheric refraction (lifts the image)
  let refractionCorrection = 0;
  if (observerState.observerOnEarth) {
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
    this.simulationTimeAtPause = Date.now(); // VERNAL_EQUINOX_2025.getTime();
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

const simStepData = { azEl: [], latLon: [] };

const earthPosWorldVec = new THREE.Vector3();
const moonPosWorldVec = new THREE.Vector3();
const embPosWorldVec = new THREE.Vector3();
const toWorldQuat = new THREE.Quaternion();

// Propagate source data at Epoch to current simulation time
// according to known rotation speeds (change rates)
function propagateOrbitalParameters(elapsedSeconds, data, changeRate) {

  const T0 = (MASTER_EPOCH.getTime() - J2000_EPOCH.getTime()) / 1000;
  const timeSinceJ2000 = T0 + elapsedSeconds;

  const julianCenturies = timeSinceJ2000 / (100 * JULIAN_YEAR);

  const current = {};
  Object.keys(data).forEach(key => {
    current[key] = data[key] + changeRate[key] * julianCenturies;

    // The mean anomaly and other parameters can grow to a very large number.
    // It's good practice to normalize it to the 0-2PI range.
    if (key !== 'A' && key !== 'EC') {
      let angle = current[key] % (2 * Math.PI);
      if (angle < 0) { angle += 2 * Math.PI; }
      current[key] = angle;
    }
  });

  return current;
}

// Animation loop function
function animate(simulation) {

  // Three.js rotation order is 'XYZ' but astronomical standard defines
  // the orientation of an orbital plane in a specific sequence:
  // 1. Rotate around the Pole
  // 2. Tilt around the Line of Nodes

  const elpasedMsec = simulation.clock.elapsed(); // in msecs
  const elapsedSeconds = elpasedMsec / 1000;

  // Propagate source data at Epoch to current simulation time
  // according to known rotation speeds (change rates)
  const currentEarth = propagateOrbitalParameters(elapsedSeconds, EARTH_J2000, EARTH_CHANGE_RATES);
  const currentMoon = propagateOrbitalParameters(elapsedSeconds, MOON_J2000, MOON_CHANGE_RATES);

  // Initial static rotations of the Sun orbital plane
  // Set the Euler rotation order to 'YXZ' to match the
  // astronomical convention for an orbital plane.
  earthOrbitalPlane.rotation.order = 'YXZ';
  earthOrbitalPlane.rotation.x = currentEarth.IN;
  earthOrbitalPlane.rotation.y = currentEarth.W_bar;

  // EMB orbits in the Earth-Sun orbital plane
  const [xe, ye, ze] = getEBMSunPosition(currentEarth);
  embPosWorldVec.set(xe, ye, ze);

  // Use a rescaled EMB-Sun distance for rendering purposes
  if (!simulation.validateMode) {
    embPosWorldVec.normalize().multiplyScalar(SUN_LIGHT_DISTANCE);
  }
  embPivot.position.copy(embPosWorldVec);

  // Earth (apsidal) orbital plane has +X axis pointing toward
  // perihelion but we want EMB to be in world frame, because
  // Earth's axial tilt orientation depends on it
  earthOrbitalPlane.updateWorldMatrix(true, false);
  earthOrbitalPlane.getWorldQuaternion(toWorldQuat);
  embPivot.quaternion.copy(toWorldQuat.invert());

  // Inclination to Ecliptic and ascending node precession
  // in backwards (retrograde) direction. Reference frame
  // is Heliocentric Ecliptic at J2000 Epoch
  moonNodesOrbitalPlane.rotation.order = 'YXZ';
  moonNodesOrbitalPlane.rotation.x = currentMoon.IN;
  moonNodesOrbitalPlane.rotation.y = currentMoon.OM;

  // Force Moon's orbital plane origin on EMB
  embPivot.getWorldPosition(embPosWorldVec);
  moonNodesOrbitalPlane.position.copy(embPosWorldVec);

  // Absidal precession in backwards (retrograde) direction
  moonAbsidalOrbitalPlane.rotation.y = currentMoon.W;

  // Moon orbit in the Moon's absidal plane
  const [xm, ym, zm] = getMoonEBMPosition(currentMoon);
  tiltedMoon.position.set(xm, ym, zm);

  // Moon's axis tilt
  // Cassini's 3rd Law  dictates that the Moon's spin axis, its orbital pole,
  // and the Ecliptic Pole are always coplanar, so pole of ecliptic (~1.54°)
  // is always *between* the 2 poles of orbit (~5.1°) and spin (~6.68°) and
  // the Moon's spin axis is "tilted back" from its orbital plane. Moon orbital
  // plane inclination to Ecliptic is done rotating around its local X axes of
  // a positive angle, so we rotate around X axis with an opposite sign.
  tiltedMoon.rotation.x = -toRadians(MOON_AXIAL_TILT);

  // Moon rotation around its Pole axis, in tidal locking with Earth
  // Moon texture map is loaded centered on 0° longitude, the Prime Meridian of
  // the Moon is aligned with its local +X axis.
  // Moon's absidal local frame +X axis direction goes from the EMB the Perigee,
  // so when the Moon is at the Perigee the 0° longitude will be on the far side
  // of the Moon. We rotate around Y axis of 180 degrees to move it facing Earth.
  const moonRotationSpeed = 2 * Math.PI / MOON_ORBITAL_PERIOD; // in tidal lock
  moon.rotation.y = Math.PI + moonRotationSpeed * elapsedSeconds;

  // Earth's axial tilt and retrograde (clockwise) precession.
  // Tilt is a negative rotation on the X-axis (Vernal Point direction). This orients
  // the North Pole towards the Winter Solstice direction, correctly setting up
  // the seasons. Precession is a slow, negative rotation on the new, tilted Y-axis.
  const timeSinceJ2000 = (MASTER_EPOCH.getTime() - J2000_EPOCH.getTime()) / 1000; // in secs
  const earthPrecessionSpeed = 2 * Math.PI / EARTH_AXIAL_PERIOD;
  tiltedEarth.rotation.x = -toRadians(EARTH_AXIAL_TILT);
  tiltedEarth.rotation.y = -earthPrecessionSpeed * (elapsedSeconds + timeSinceJ2000);

  // Earth daily rotation
  // Earth texture map is loaded centered on 0° longitude, the Prime Meridian of
  // the Earth is aligned with its local +X axis that is also the world +X axis
  // and points from Earth toward Vernal Point. Greenwich Apparent Sidereal
  // Time (GAST) is the angle between the Vernal Point and the Prime Meridian at
  // master epoch time.
  const earthRotationSpeed = 2 * Math.PI / SIDERAL_DAY;
  const gast = EARTH_DATA.L_Ap_Sid_Time ;
  earth.rotation.y = gast + earthRotationSpeed * elapsedSeconds;

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

  if (simulation.validateMode) {
    // Force hierarchy recalculation and recompute real positions
    scene.updateWorldMatrix(true, true);

    const m = tiltedMoon.getWorldPosition(moonPosWorldVec);
    const e = tiltedEarth.getWorldPosition(earthPosWorldVec);
    const emb = embPivot.getWorldPosition(embPosWorldVec);

    const em = m.clone().sub(emb); // EMB -> Moon (in world frame)
    offsetToLocal(em, embPivot); // EMB -> Moon (in EMB frame)

    const moonV  = convertToJPL(em);
    const earthV = convertToJPL(e);
    const timeV = elpasedMsec + simulation.clock.masterEpochTime;
    return {
      time: timeV,         // msec
      sunPosition: earthV, // km
      moonPosition: moonV  // km
    };
  }

  // Use a rescaled EMB-Sun distance for rendering purposes
  //embPivot.position.normalize().multiplyScalar(SUN_LIGHT_DISTANCE);
  scene.updateWorldMatrix(true, true);

  // Set Sun light to look at Earth
  tiltedEarth.getWorldPosition(earthPosWorldVec);
  sunLight.target.position.copy(earthPosWorldVec);

  // Update the end point of the line to match the sun's new position
  if (DEBUG) {
    set_sunline_length(earthPosWorldVec);
    //embPivot.getWorldPosition(embPosWorldVec);
    //scene_axes.position.copy(embPosWorldVec);
  }

  if (observerState.object) {
    // In dry run mode sample the proper value and return
    if (simulation.dryRunFunction) {
        return simulation.dryRunFunction();
    }
    // Update celestial events dependent on observer position
    const timeForward = (simulation.clock.speed() > 0);
    const sunDistance = earthPosWorldVec.length();
    simulation.eventManager.update(timeForward, sunDistance, tiltedMoon, sunLight);

    // If we are in observer view pass elevation and azimuth
    const view = views.get(observerState.viewIndex);
    const isAct = (view == views.getActive());
    const { azEl, latLon } = isAct ? view.getGeoData(observerState.marker) : { azEl: [], latLon: []};
    Object.assign(simStepData, { azEl: azEl, latLon: latLon });
  }

  // If is first frame init main view and other stuff that requires
  // updated orbital parameters and Earth position
  if (views.getActive() === null) {
    views.init(earthPosWorldVec);
    init_moon_path(currentMoon);
  }

  // Update the active controls and render with the active camera
  const activeCamera = views.update();
  renderer.render(scene, activeCamera);

  return simStepData;
}

// Functions used in dry run to sample the requested sim output value
function sunHeight() {
  // Calculate the equivalent radius using the captured real distance
  const sunDistance = earthPosWorldVec.length();
  const equivalentSunRadius = toUnits(SUN_RADIUS_KM) * (SUN_LIGHT_DISTANCE / sunDistance);
  return heightAboveHorizon(sunLight, equivalentSunRadius);
}

function moonHeight() {
  return heightAboveHorizon(tiltedMoon, toUnits(MOON_RADIUS_KM));
}

// ============================================================================
// EXPORT CLASS
// ============================================================================

class Simulation {
  constructor(validateMode) {
    // Don't render and return Moon / Earth positions
    this.validateMode = validateMode;
    this.dryRunFunction = null;
    this.dryRunFunctions = [sunHeight, sunHeight, moonHeight, moonHeight];

    // Our simulation clock
    this.clock = new SimClock(MASTER_EPOCH);

    // Create a new event manager and expose its 'on' method
    this.eventManager = new CelestialEventManager(this, heightAboveHorizon);
    this.on = this.eventManager.on.bind(this.eventManager);

    this.update = () => animate(this);
    this.getTime = this.clock.getTime.bind(this.clock);
    this.setDate = this.clock.setDate.bind(this.clock);
    this.speed = this.clock.speed.bind(this.clock);
    this.setSpeed = this.clock.setSpeed.bind(this.clock);
    this.togglePause = this.clock.togglePause.bind(this.clock);
    this.findNextEvent = this.eventManager.findNextEvent.bind(this.eventManager);
    this.setActiveView = views.setActive.bind(views);
    this.disposeView = views.dispose.bind(views);
    this.lockToOrbit = lockToOrbit;
    this.enterObserverMode = enterObserverMode;
    this.placeObserverAt = placeObserverAt;
    this.exitObserverMode = exitObserverMode;
    this.pickObject = pickObject;
  }
  setDryRunFunction(idx) {
    this.dryRunFunction = (idx === null ? null : this.dryRunFunctions[idx]);
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
    views.setDefault();
  }
  isObserverView() {
    return views.getActive() === views.get(observerState.viewIndex);
  }
  isOrbitLocked(object) {
    const lockedObjects = views.getOrbitLockedObjects();
    return lockedObjects.includes(object);
  }
  unlockCamera() {
    const view = views.getActive();
    const marker = view.cameraLock.target;
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
