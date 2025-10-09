'use strict';

// Physical properties (in kilometers, converted to scene units)
const EARTH_RADIUS_KM = 6371; // Conventional radius for sphere approx.
const MOON_RADIUS_KM = 1737.53;
const SUN_RADIUS_KM = 695700;
const MOON_DISTANCE_KM = 384400;    // average

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

  // We always return in [x1, x2] order. This ensures (start < end) == (x1 < x2)

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

export class CelestialEventManager {
  constructor(simulation,visibilityFunction) {
    this.listeners = new Map();
    this.sim = simulation;
    this.visibilityFunction = visibilityFunction;
    this.state = {};
  }
  on(eventType, callback) {
    if (!this.listeners.has(eventType)) {
        this.listeners.set(eventType, []);
    }
    this.listeners.get(eventType).push(callback);
  }
  emit(eventType, data) {
    const listeners = this.listeners.get(eventType) || [];
    listeners.forEach(callback => callback(data));
  }

  checkVisibility(object, name, radius, timeForward, item) {

    const isVisible = this.visibilityFunction(object, radius) > 0;

    const state = this.state;
    if (!state.hasOwnProperty(item)) { // First frame
        state[item] = isVisible;
        return isVisible;
    };

    if (isVisible !== state[item]) {
        const atRise = (isVisible === timeForward);
        if (atRise) {
            this.emit(CELESTIAL_EVENTS.RISE, name);
        } else {
            this.emit(CELESTIAL_EVENTS.SET, name);
        }
        state[item] = isVisible;
    }
    return isVisible;
  }

  update(timeForward, sunDistance, tiltedMoon, sunLight) {
      // Calculate the equivalent radius using the captured real distance
      const equivalentSunRadius = toUnits(SUN_RADIUS_KM) * (SUN_LIGHT_DISTANCE / sunDistance);

      // Check and emit celestial events
      this.checkVisibility(
          tiltedMoon, "Moon", toUnits(MOON_RADIUS_KM), timeForward, "moonVisible"
      );
      this.checkVisibility(
          sunLight, "Sun", equivalentSunRadius, timeForward, "sunVisible"
      );
  }

  doSimStepAt(currentTimeMs) {
    const currentDate = new Date(currentTimeMs);
    this.sim.setDate(currentDate); // Set sim date, UI will read from it
    return this.sim.update();
  }

  findBracket(start, forward, first_sign, tolerance) {

    const NINE_HOURS = 7 * 3600 * 1000; // 9 hours interval (in msec)
    const step = (forward ? NINE_HOURS : -NINE_HOURS);
    let current = start + step / 18; // Start some minutes before/after now
    let currentValue = this.doSimStepAt(current);
    let next = current + step;

    if (Math.abs(currentValue) < tolerance) {
      current -= (current - start) / 2; // jump back a bit, but still after start
      currentValue = this.doSimStepAt(current);
    }

    while (true) {
      const nextValue = this.doSimStepAt(next);

      if (Math.abs(nextValue) < tolerance) {
        next += step / 18; // Jump a bit forward
        continue;
      }

      if (currentValue * nextValue < 0 &&  // Different sign with next
          currentValue * first_sign > 0) { // Same sign with first_sign
          return [current, next, currentValue, nextValue];
      }

      current = next;
      next += step;
      currentValue = nextValue;
    }
  }

  findNextEvent(event, forward) {
    // Set function that will be called in Simulation to pick
    // and return the specific value according to event type
    const idx = CELESTIAL_EVENTS.names.indexOf(event);
    this.sim.setDryRunFunction(idx);

    const isRise = CELESTIAL_EVENTS.types[idx] === "rise";
    const first_sign = (isRise === forward ? -1 : 1);

    // Find initial [star, end] range
    const tolerance = 0.01;
    const now = this.sim.getTime(); // in msec
    const [startTime, endTime, startValue, endValue] = this.findBracket(now, forward, first_sign, tolerance);

    // Solve the position, zero is within returned range with x2
    // the the timestamp always _after_ the event
    const [x1, x2] = linearSolver(startTime, endTime, startValue, endValue, this.doSimStepAt.bind(this), tolerance);

    this.sim.setDryRunFunction(null);
    this.state = {}; // Reset stale states
    return forward ? x2 : x1;
  }
};
