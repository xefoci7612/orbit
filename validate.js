'use strict';

// JPL Validation Script

// JPL settings are:
// Target Body: Sun [Sol]
// Coordinate Center: 0°E, 0°N, 0 k
// Time Specification: startDate to endDate, stepInDays
// Table Settings: defaults

import * as THREE from 'three';

// Runs the simulation over a date range and triggers a CSV download
export async function runValidation(sim) {
  // All configuration is in this single object for easy access
  const config = {
    startDate: new Date('2025-09-10T09:34:03Z').toISOString(),
    endDate: new Date('2025-12-10T09:34:03Z').toISOString(),
    stepInDays: 1,
    // How many simulation steps to process before pausing to keep the browser responsive
    chunkSize: 100,
  };

  console.log('Starting validation...');
  sim.validateMode = true;

  try {
    // Generate the CSV content asynchronously
    const csvContent = await generateJplCsvContent(sim, config);

    // Create a filename and trigger the download
    const fileName = `validation_${config.startDate.slice(0, 10)}_${config.endDate.slice(0, 10)}.csv`;
    downloadFile(csvContent, fileName, 'text/csv');

    console.log('Validation complete. File download initiated.');

  } catch (error) {
    console.error('An error occurred during validation:', error);
  } finally {
    sim.validateMode = false;
  }
}

// Converts a position vector from the simulation's internal units and coordinate
// system to the JPL Horizons Ecliptic J2000 standard (in kilometers).
function convertToJPL(simPositionInKm) {
  const [sim_x_km, sim_y_km, sim_z_km] = simPositionInKm; // sim data is in km

  // This transformation maps the simulation's world frame to the JPL Ecliptic J2000 frame
  const jpl_x =   sim_x_km; // Matches directly
  const jpl_y = - sim_z_km; // Negate the Z component
  const jpl_z =   sim_y_km; // Use the Y component (which is 0 for the Sun)

  return [jpl_x, jpl_y, jpl_z];
}

// Iterates through simulation steps and generates JPL-formatted CSV data
async function generateJplCsvContent(sim, config) {
  const startMs = Date.parse(config.startDate);
  const endMs = Date.parse(config.endDate);
  const stepMs = config.stepInDays * 24 * 3600 * 1000;

  // Updated header to reflect the new coordinate system
  const csvRows = ['ISO,jpl_sun_x,jpl_sun_y,jpl_sun_z,jpl_moon_x,jpl_moon_y,jpl_moon_z'];
  let stepsProcessed = 0;

  for (let currentTimeMs = startMs; currentTimeMs <= endMs; currentTimeMs += stepMs) {
    const currentDate = new Date(currentTimeMs);
    sim.setDate(currentDate);

    // Get the raw data from the simulation
    const rawData = sim.update();

    // Convert positions to the JPL coordinate system
    const jplSunPosition = convertToJPL(rawData.sunPosition, true);
    const jplMoonPosition = convertToJPL(rawData.moonPosition);

    // Format the CONVERTED data into a CSV row
    const row = [
      currentDate.toISOString(),
      ...jplSunPosition,
      ...jplMoonPosition
    ].join(',');
    csvRows.push(row);

    stepsProcessed++;
    if (stepsProcessed % config.chunkSize === 0) {
      await new Promise(resolve => setTimeout(resolve, 0));
    }
  }

  return csvRows.join('\n');
}

// Trigger a browser download for the given content
function downloadFile(content, fileName, mimeType) {
  const blob = new Blob([content], { type: mimeType });
  const a = document.createElement('a');

  a.href = URL.createObjectURL(blob);
  a.download = fileName;
  a.click();

  // Clean up by revoking the object URL to free up memory
  URL.revokeObjectURL(a.href);
}


// Helper debug functions

// Physical properties (in kilometers, converted to scene units)
const EARTH_RADIUS_KM = 6371; // Conventional radius for sphere approx.
const MOON_RADIUS_KM = 1737.53;
const SUN_RADIUS_KM = 695700;
const MOON_DISTANCE_KM = 384400;    // average

// Scale conversions
const KM_PER_UNIT = EARTH_RADIUS_KM / 50; // Earth radius to units
const toUnits = km => km / KM_PER_UNIT;

function addAxesHelper(object, size) {
    const axesHelper = new THREE.AxesHelper(size);
    axesHelper.setColors(
      new THREE.Color(0xff0000), // Red
      new THREE.Color(0x00ff00), // Green
      new THREE.Color(0x0000ff)  // Blue
    );
    object.add(axesHelper);
    return axesHelper;
}

function createMeridianLine(radius, longitudeRad = 0, segments = 64, color = 0xffff00) {
  const points = [];

  // Loop from North Pole to South Pole
  for (let i = 0; i <= segments; i++) {
    // Interpolate the latitude from +90 degrees to -90 degrees
    const latRad = Math.PI / 2 - (i / segments) * Math.PI;

    // Calculate the 3D position using spherical coordinates (with Y as up)
    const x = radius * Math.cos(latRad) * Math.cos(longitudeRad);
    const y = radius * Math.sin(latRad);
    const z = radius * Math.cos(latRad) * Math.sin(longitudeRad);

    points.push(new THREE.Vector3(x, y, z));
  }

  const geometry = new THREE.BufferGeometry().setFromPoints(points);
  const material = new THREE.LineBasicMaterial({ color: color, fog: false }); // fog:false keeps it visible from far away
  const line = new THREE.Line(geometry, material);

  return line;
}

// Create a bright yellow line to represent the Prime Meridian (0° longitude).
// We make the radius slightly larger than the Earth's to prevent Z-fighting
// (where the line flickers because it's at the same depth as the texture).
const primeMeridianLine = createMeridianLine(toUnits(EARTH_RADIUS_KM) * 1.01);

// Create a red material for the line
const sunLineMaterial = new THREE.LineBasicMaterial({ color: 0xffff00 });

// Define the start and end points for the line. Will be updated by animation
const sunLinePoints = [];
sunLinePoints.push(new THREE.Vector3(0, 0, 0)); // Start point
sunLinePoints.push(new THREE.Vector3(0, 0, 0)); // End point

const sunLineGeometry = new THREE.BufferGeometry().setFromPoints(sunLinePoints);

// Create the line object
const sunLine = new THREE.Line(sunLineGeometry, sunLineMaterial);

export function init_debug(scene, tiltedEarth, earth) {
  earth.add(primeMeridianLine);
  scene.add(sunLine); // in world coordinates
  addAxesHelper(tiltedEarth, toUnits(2 * EARTH_RADIUS_KM));
}

// Update the sunline to match the earth new position
export function set_sunline_length(earthPosWorldVec) {
  const s = earthPosWorldVec.clone();
  const d = toUnits(MOON_DISTANCE_KM / 2);
  const offset = s.clone().normalize().multiplyScalar(d).negate();
  const e = s.clone().add(offset);
  const positions = sunLine.geometry.attributes.position.array;
  s.toArray(positions, 0);
  e.toArray(positions, 3);

  // Tell Three.js that the position attribute needs to be updated on the GPU
  sunLine.geometry.attributes.position.needsUpdate = true;
}
