'use strict';

// Constants
const GM = 398600.4418; // Earth gravitational parameter (km^3/s^2)
const SECONDS_PER_DAY = 86400;

// Scale conversions
const EARTH_RADIUS_KM  = 6371;
const KM_PER_UNIT = EARTH_RADIUS_KM / 50; // Earth radius to units
const toUnits = km => km / KM_PER_UNIT;
const toRadians = degrees => degrees * Math.PI / 180;

const SAT_GROUPS = [
  'Stations',
  'Geo',
  'Intelsat',
  'Eutelsat',
  'Geodetic',
];

const SAT_DATA_READY = 'satellites_ready';
const SAT_CHANGED = 'satellite_changed';

/* Sample satellite data format from Celestrak (JSON):

    OBJECT_NAME	"EUTELSAT 7A"
    OBJECT_ID	"2004-008A"
    EPOCH	"2025-11-04T13:54:33.454368"
    MEAN_MOTION	1.00273851
    ECCENTRICITY	0.0001438
    INCLINATION	5.2575
    RA_OF_ASC_NODE	74.6779
    ARG_OF_PERICENTER	184.6037
    MEAN_ANOMALY	214.2663
    EPHEMERIS_TYPE	0
    CLASSIFICATION_TYPE	"U"
    NORAD_CAT_ID	28187
    ELEMENT_SET_NO	999
    REV_AT_EPOCH	7910
    BSTAR	0
    MEAN_MOTION_DOT	1.17e-6JS:0.00000117
    MEAN_MOTION_DDOT	0
*/

async function fetchSatellites(satGroup, referenceTime) {

  const encSatGroup = encodeURIComponent(satGroup)
  const url = `https://celestrak.org/NORAD/elements/gp.php?GROUP=${encSatGroup}&FORMAT=json`;

  try {
    // Fetch satellite data
    const response = await fetch(url);
    if (!response.ok) throw new Error(`HTTP error: ${response.status}`);

    const satellitesData = await response.json();

    // Process each satellite
    return satellitesData.map(sat => {

      // Calculate semi-major axis from mean motion (revolutions/day) and
      // Kepler's third law
      const meanMotionRadPerSec = (sat.MEAN_MOTION * 2 * Math.PI) / SECONDS_PER_DAY;
      const semiMajorAxisKm = Math.pow(GM / (meanMotionRadPerSec ** 2), 1/3);
      const satEpoch = new Date(sat.EPOCH);

      // Create orbital elements object
      return {
        NAME: sat.OBJECT_NAME,
        A: toUnits(semiMajorAxisKm),          // Semi-major axis in scene units
        EC: sat.ECCENTRICITY,                 // Eccentricity
        IN: toRadians(sat.INCLINATION),      // Inclination (rad)
        MA: toRadians(sat.MEAN_ANOMALY),     // Mean Anomaly (rad)
        W: toRadians(sat.ARG_OF_PERICENTER), // Argument of Perigee (rad)
        OM: toRadians(sat.RA_OF_ASC_NODE),   // RAAN (rad)
        EPOCH: sat.EPOCH,                    // Epoch as Date object
        TIME: (satEpoch.getTime() - referenceTime) / 1000, // in secs
        NORAD_ID: sat.NORAD_CAT_ID,          // For reference
      };
    });

  } catch (error) {
    console.error('Satellite data fetch failed:', error);
    return []; // Return empty array on failure
  }
}

/**
 * Satellite data manager with event system
 */

class Satellites {
  constructor(referenceDate) {
    this.referenceTime = referenceDate.getTime();
    this.listeners = new Map(); // Event listeners
    this.satellitesData = new Map(); // Cached satellite data
    this.isFetching = false;
    this.activeGroup = null;
    this.activeSat = null;
  }

  getActiveSat() {
    return this.activeSat;
  }

  getActiveGrop() {
    return this.activeGroup;
  }

  setActiveSat(satName) {
    const group = this.activeGroup;
    const sats = group ? this.satellitesData.get(group) : [];
    const sat = sats.find(sat => sat.NAME === satName);
    if (!sat) {
      throw new Error("Satellite does not belong to the active group");
    }
    if (!this.activeSat || this.activeSat.NAME !== satName) {
      this.activeSat = sat;
      this.emit(SAT_CHANGED, sat);
    }
  }

  setActiveGroup(satGroup) {
    if (!this.satellitesData.has(satGroup)) {
      this.loadSatellites(satGroup);
      return;
    }
    this.activeGroup = satGroup;
    const satellites = this.satellitesData.get(satGroup);
    const name = satellites[0].NAME;
    this.setActiveSat(name); // Set first satellite as active
    this.emit(SAT_DATA_READY);
  }

  get(satGroup) {
    return this.satellitesData.get(satGroup) || [];
  }

  on(eventType, callback) {
    if (!this.listeners.has(eventType)) {
      this.listeners.set(eventType, []);
    }
    this.listeners.get(eventType).push(callback);
  }

  emit(eventType, data) {
    const callbacks = this.listeners.get(eventType) || [];
    callbacks.forEach(callback => callback(data));
  }

  onDataReady(satellites, satGroup) {
    this.isFetching = null;
    this.satellitesData.set(satGroup, [...satellites]);
    if (satellites.length > 0) {
      this.setActiveGroup(satGroup);
    }
  }

  loadSatellites(satGroup) {
    if (!this.isFetching && !this.satellitesData.has(satGroup)) {
      this.isFetching = true;
      fetchSatellites(satGroup, this.referenceTime).then(data => {
        this.onDataReady(data, satGroup);
      })
      .catch(error => {
        console.error('Failed to load satellites:', error);
        this.isFetching = false; // allow retry if needed
      });
    }
  }
};

export { Satellites, SAT_DATA_READY, SAT_CHANGED, SAT_GROUPS };
