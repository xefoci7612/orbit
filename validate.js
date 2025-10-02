'use strict';

// Set this to false to run the interactive UI instead of the validation script.
export const Validate = true;

// Runs the simulation over a date range and triggers a CSV download
export async function runValidation(sim) {
  // All configuration is in this single object for easy access
  const config = {
    startDate: '2025-09-10T00:00:00Z',
    endDate: '2025-10-10T00:00:00Z',
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
