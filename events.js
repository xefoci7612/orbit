'use strict';

// Physical properties (in kilometers, converted to scene units)
const EARTH_RADIUS_KM = 6371; // Conventional radius for sphere approx.
const MOON_RADIUS_KM = 1737.53;
const SUN_RADIUS_KM = 695700;
const MOON_DISTANCE_KM = 384400; // average

// Scale conversions
const KM_PER_UNIT = EARTH_RADIUS_KM / 50; // Earth radius to units
const toUnits = km => km / KM_PER_UNIT;

// Rescaled Sun distance for placing directional Sun light source
const SUN_LIGHT_DISTANCE = toUnits(10 * MOON_DISTANCE_KM);

// Event types
export const CELESTIAL_EVENTS = {
    RISE: 'rise',
    SET: 'set',

    names: ["sunrise", "sunset", "moonrise", "moonset"],
    types: ["rise", "set", "rise", "set"],
};


// Narrows down a time bracket using the regula falsi method and returns a [x1, x2]
// high-precision bracket. The order is preserved, if start < end, then t1 < t2.
function linearSolver(start, end, fStart, fEnd, fun, tolerance, maxIterations = 10) {
  let x1 = start; // times in msec
  let x2 = end;
  let f1 = fStart;
  let f2 = fEnd;

  // Verify that the function has opposite signs at the endpoints
  if (f1 * f2 > 0) {
      throw new Error("Function must have opposite signs at start and end points");
  }

  for (let i = 0; i < maxIterations; i++) {

      // Check if interval is smaller than tolerance
      if (Math.abs(x2 - x1) < tolerance) {
          return [x1, x2];
      }

      // Calculate the next guess using linear interpolation (regula falsi method)
      const xNext = (x1 * f2 - x2 * f1) / (f2 - f1);
      const fNext = fun(xNext);

      // Check if we've found the root within tolerance
      if (Math.abs(fNext) < tolerance) {
          return [x1, x2];
      }

      // Update the interval based on the sign of fNext
      if (f1 * fNext < 0) {
          // Root is between x1 and xNext
          x2 = xNext;
          f2 = fNext;
      } else {
          // Root is between xNext and x2
          x1 = xNext;
          f1 = fNext;
      }

      if (i + 1 === maxIterations)
        console.log("Max iterations reached in linearSolver");
  }

  // Return the best estimate after max iterations
  return [x1, x2];
}

// Main class to handle celestial events functionality
export class CelestialEventManager {
  constructor(simulation, fun) {
    this.listeners = new Map();
    this.sim = simulation;
    this.heightAboveHorizon = fun;
    this.prevState = {};
  }

  // Register a callback for a given event type
  on(eventType, callback) {
    if (!this.listeners.has(eventType)) {
        this.listeners.set(eventType, []);
    }
    this.listeners.get(eventType).push(callback);
  }

  // Call the registered callbacks on event emit
  emit(eventType, data) {
    const listeners = this.listeners.get(eventType) || [];
    listeners.forEach(callback => callback(data));
  }

  // Check for visibility change, update internal state,
  // and return the event to emit
  checkNewEvent(object, radius, item, timeForward) {
    // Call the external function that performs the actual
    // visibility check in context
    const isVisible = this.heightAboveHorizon(object, radius) > 0;
    const atRise = (isVisible === timeForward);
    const event = atRise ? CELESTIAL_EVENTS.RISE : CELESTIAL_EVENTS.SET;
    const firstFrame = !this.prevState.hasOwnProperty(item);
    const visibilityChanged = (!firstFrame && isVisible !== this.prevState[item]);
    this.prevState[item] = isVisible;
    return visibilityChanged ? event : null;
  }

  // Emit event if visibility changed since last frame
  update(timeForward, sunDistance, tiltedMoon, sunLight) {
    // Calculate the equivalent radius using the captured real distance
    const sunRadius = toUnits(SUN_RADIUS_KM) * (SUN_LIGHT_DISTANCE / sunDistance);
    const moonRadius = toUnits(MOON_RADIUS_KM);
    const eventM = this.checkNewEvent(tiltedMoon, moonRadius, "moonVisible", timeForward);
    if (eventM !== null) {
      this.emit(eventM, "Moon");
    }
    const eventS = this.checkNewEvent(sunLight, sunRadius, "sunVisible", timeForward);
    if (eventS !== null) {
      this.emit(eventS, "Sun");
    }
  }

  // Perform a simulation dryrun step after setting simulation clock to given time
  doSimStepAt(currentTimeMs) {
    const currentDate = new Date(currentTimeMs);
    this.sim.setDate(currentDate); // Update sim clock, UI will read from it
    return this.sim.update();
  }

  // Perform a fixed step search for a suitable range [x1, x2] that contains the
  // correct type of event (rise/set). Range will be the input for linear solver,
  // so we must ensure signs of values at f(x1) and f(x2) are different
  findBracket(start, step, first_sign, tolerance) {
    let x1 = start;
    let x2 = x1 + step;
    let f1 = this.doSimStepAt(x1);

    // If we are on a zero, jump a bit forward. Our job
    // is to find a suitable range, not a zero
    if (Math.abs(f1) < tolerance) {
      x1 += step / 18;
      f1 = this.doSimStepAt(x1);
    }

    // For now leave the search loop unguarded
    while (true) {

      const f2 = this.doSimStepAt(x2);

      // If we are on a zero, jump a bit forward and resample
      if (Math.abs(f2) < tolerance) {
        x2 += step / 18;
        continue;
      }

      // Ensure f1 sign is correct to properly differentiate
      // rise (negative->positive) and set (positive->negative)
      if (f1 * f2 < 0 &&  // Different signs
          f1 * first_sign > 0) { // Same sign with first_sign
          return [x1, x2, f1, f2];
      }

      x1 = x2;
      x2 += step;
      f1 = f2;
    }
  }

  // Find the next event in time direction starting from now.
  // Select a specific external function and use a linear solver
  // to find the zero, actually a very small range around zero
  findNextEvent(event, forward) {
    // Set the external function that will be called in Simulation
    // loop to pick the specific value according to event type
    const idx = CELESTIAL_EVENTS.names.indexOf(event);
    this.sim.setDryRunFunction(idx);

    // Find a zero is not enough, we have to ensure the 'direction'
    // of the function, so to disambiguate between rise and set events
    const isRise = CELESTIAL_EVENTS.types[idx] === "rise";
    const first_sign = (isRise === forward ? -1 : 1); // Sign of range's first endpoint

    // Chose a proper time step for findBracket.
    // Step should be large enough for efficency but also smaller than the minimum
    // interval between two consecutive (signed) zeros of the function.
    const NINE_HOURS = 9 * 3600 * 1000; // in msec
    const step = (forward ? NINE_HOURS : -NINE_HOURS);

    // Find initial [x1, x2] range
    const tolerance = 0.01;
    const start = this.sim.getTime() + step / 18; // Start some minutes before/after now
    const [x1, x2, f1, f2] = this.findBracket(start, step, first_sign, tolerance);

    // Find the zero of the function, that must be monotonic in [x1, x2], and return a
    // high-precision range [z1, z2] satisfying the condition (x1 < x2) == (z1 < z2)
    const [z1, z2] = linearSolver(x1, x2, f1, f2, this.doSimStepAt.bind(this), tolerance);

    this.sim.setDryRunFunction(null);
    this.prevState = {}; // Reset stale states
    return forward ? z2 : z1; // Pick the timestamp _after_ the event
  }
};
