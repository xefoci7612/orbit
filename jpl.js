'use strict';

// File names
const JPL_EMB_SSB_FILE = "jpl_emb_sun.txt"

// Scale conversions
const EARTH_RADIUS_KM = 6371; // Conventional radius for sphere approx.
const KM_PER_UNIT = EARTH_RADIUS_KM / 50; // Earth radius to units
const toUnits = km => km / KM_PER_UNIT;
const toRadians = degrees => degrees * Math.PI / 180;

// The Julian Date of the Unix Epoch (a constant)
const JD_OF_UNIX_EPOCH = 2440587.5;

// Convert JD to milliseconds since the Unix Epoch
const toMsec = jd => (jd - JD_OF_UNIX_EPOCH) * 86400 * 1000;

// Data records
let emb_ssb_data = null;

async function loadAndParseJplRecords(filename) {

    let response = null;
    try {
      response = await fetch(filename);
      if (!response.ok) {
        throw new Error(`Failed to load file: ${response.status} ${response.statusText}`);
      }
    } catch (error) {
      console.error('Error loading or parsing JPL data:', error);
      throw error;
    }

    const text = await response.text();
    const lines = text.split('\n');
    const records = [];

    // Get only the lines between the start and end markers
    const dataLines = lines.slice(
      lines.indexOf('$$SOE') + 1,
      lines.indexOf('$$EOE')
    );
    /*
     Format of each record is:

        2460928.898645833 = A.D. 2025-Sep-10 09:34:03.0000 TDB
         EC= 2.348399511484001E-02 QR= 1.451551996030375E+08 IN= 9.565607274919776E-03
         OM= 2.345912713353362E+02 W = 2.311018377459241E+02 Tp=  2461045.467873353045
         N = 1.152491697491665E-05 MA= 2.439258621986873E+02 TA= 2.415401694487883E+02
         A = 1.486460015779342E+08 AD= 1.521368035528309E+08 PR= 3.123666754246648E+07
    */
    const patterns = [
      /=\s*A\.D\.\s*([-\w\s:]+)\.\d+\s*TDB/, // /^([\d.]+)/, // Julian Date
      /EC=\s*([\d.E+-]+)\s*QR=\s*(?:[\d.E+-]+)\s*IN=\s*([\d.E+-]+)/,
      /OM=\s*([\d.E+-]+)\s*W =\s*([\d.E+-]+)\s*Tp=\s*(?:[\d.E+-]+)/,
      /MA=\s*([\d.E+-]+)\s*TA=\s*(?:[\d.E+-]+)/,
      /A =\s*([\d.E+-]+)\s*AD=\s*(?:[\d.E+-]+)\s*PR=\s*([\d.E+-]+)/
    ];

    // Process the data in chunks of 5 lines
    for (let i = 0; i < dataLines.length; i += 5) {
      const recordLines = dataLines.slice(i, i + 5);
      const allMatches = [];
      for (let j = 0; j < recordLines.length; j += 1) {
        const line = recordLines[j];
        const pattern = patterns[j];
        allMatches.push(line.match(pattern));
      }
      const item = {
        date: toMsec(parseFloat(allMatches[0][1])),
        EC: parseFloat(allMatches[1][1]),
        IN: toRadians(parseFloat(allMatches[1][2])),
        OM: toRadians(parseFloat(allMatches[2][1])),
        W: toRadians(parseFloat(allMatches[2][2])),
        MA: toRadians(parseFloat(allMatches[3][1])),
        A: toUnits(parseFloat(allMatches[4][1])),
        PR: parseFloat(allMatches[4][2]),
      };
      records.push(item);
    }

    console.log(`Successfully loaded ${records.length} JPL records.`);

    return records;
}

export async function init_jpl_data() {

 emb_ssb_data = await loadAndParseJplRecords(JPL_EMB_SSB_FILE);
}

// Find the index of the first record whose date is after the requested time
function find_closest_pair(time) {

  const upperIndex = emb_ssb_data.findIndex(record => record.date > time);

  // Handle edge cases
  if (upperIndex === -1) { // after last record
    return { before: emb_ssb_data[emb_ssb_data.length - 1], after: null };
  } else if (upperIndex === 0) { // before first record
    return { before: null, after: emb_ssb_data[0] };
  }

  return {
    before: emb_ssb_data[upperIndex - 1],
    after: emb_ssb_data[upperIndex]
  };
}

function interpolate(before, after, time) {

  // Guard against null inputs or identical timestamps
  if (!before || !after) return before || after || null;

  // Define which fields are angles (in radians) and need special handling
  const angularKeys = new Set(['IN', 'OM', 'W', 'MA']);
  const TWO_PI = 2 * Math.PI;

  // Calculate the interpolation factor (t)
  const totalTimeSpan = after.date - before.date;
  const fraction = (time - before.date) / totalTimeSpan;

  const interpolatedRecord = { date: time };

  for (const key in before) {

    if (key === 'date') {
      continue;
    }

    const beforeVal = before[key];
    const afterVal = after[key];

    // Angles require handling of 0/2PI boundary
    if (angularKeys.has(key)) {
      // If the difference is more than 180° the shortest path is the other way
      // around the circle, so subtract 360° from the difference
      let diff = afterVal - beforeVal;
      if (diff > Math.PI) {
        diff -= TWO_PI;
      } else if (diff < -Math.PI) {
        diff += TWO_PI;
      }

      const value = beforeVal + fraction * diff;

      // Normalize the result to be within [0, 2PI)
      interpolatedRecord[key] = (value % TWO_PI + TWO_PI) % TWO_PI;

    } else {
      interpolatedRecord[key] = beforeVal + fraction * (afterVal - beforeVal);
    }
  }

  return interpolatedRecord;
}

// Get the ephemeris data for a specific time, interpolating if necessary
export function getEMB_data(time) {
  const pair = find_closest_pair(time);
  return interpolate(pair.before, pair.after, time);
}
