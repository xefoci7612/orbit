'use strict';

// Define events and associated type/function
// period is in days, direction disambiguates between rise and set events
const eventDefinitions = [
  { name: 'sunrise',    type: 'sunRiseSet', needObserver: true, period: 1, ascending: true },
  { name: 'sunset',     type: 'sunRiseSet', needObserver: true, period: 1, ascending: false},
  { name: 'moonrise',   type: 'moonRiseSet',needObserver: true, period: 1, ascending: true },
  { name: 'moonset',    type: 'moonRiseSet',needObserver: true, period: 1, ascending: false},
  { name: 'suntransit', type: 'sunTransit', needObserver: true, period: 1, ascending: true },
  { name: 'moonnew',    type: 'newMoon',    needObserver:false, period:27, ascending: true },
];

// These names is how UI html refers to a specific event
export const CELESTIAL_EVENTS = {
  event: "celestial_event",
  names: eventDefinitions.map(def => def.name),
};

// Main class to handle celestial events functionality
export class CelestialEventManager {
  constructor(simulation, eventFuns) {
    this.listeners = new Map();
    this.sim = simulation;
    this.prevState = {};

    // The external in-context event detection functions: functions must be
    // in same order as the types in eventDefinitions!
    const keys = [...new Set(eventDefinitions.map(def => def.type))];
    this.eventFunctions = Object.fromEntries(keys.map((k, i) => [k, eventFuns[i]]));

    this.eventChecks = [
      { target:  'sun', fun: this.checkRiseSetEvent.bind(this), type: 'sunRiseSet' },
      { target: 'moon', fun: this.checkRiseSetEvent.bind(this), type: 'moonRiseSet'},
      { target:  'sun', fun: this.checkTransitEvent.bind(this), type: 'sunTransit' },
      { target: 'moon', fun: this.checkTransitEvent.bind(this), type: 'newMoon'    },
    ];
  }

  reset() {
    this.prevState = {}; // stale events after observer changes
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
  checkRiseSetEvent(type, timeForward) {
    // Call the external function that performs the actual
    // visibility check in context
    const fun = this.eventFunctions[type];
    const isVisible = fun() > 0;
    const atRise = (isVisible === timeForward);
    const event = atRise ? 'rise' : 'set';
    const firstFrame = !this.prevState.hasOwnProperty(type);
    const visibilityChanged = (!firstFrame && isVisible !== this.prevState[type]);
    this.prevState[type] = isVisible;
    return visibilityChanged ? event : null;
  }

  // Check for transit event
  // We should be careful because signedDistance is a tri-state [-1, 0, 1]
  // in particular we risk double fire if signedDistance == 0, so check
  // only for a strictly positive afterMeridian.
  checkTransitEvent(type, timeForward) {
    const fun = this.eventFunctions[type];
    const signedDistance = fun();
    const hasFullyCrossed = (signedDistance > 0) === timeForward;
    const event = (type === 'newMoon' ? 'new' : 'transit');
    const firstFrame = !this.prevState.hasOwnProperty(type);
    const signChanged = (!firstFrame && signedDistance * this.prevState[type] <= 0);
    this.prevState[type] = signedDistance;
    return signChanged && hasFullyCrossed ? event : null;
  }

  // Emit new events since last frame
  update(timeForward, isObserver) {
    for (const check of this.eventChecks) {
      const item = eventDefinitions.find(item => item.type === check.type);
      if (!isObserver && item.needObserver)
        continue;
      const eventType = check.fun(check.type, timeForward);
      if (eventType !== null) {
        this.emit(CELESTIAL_EVENTS.event, check.target + eventType);
      }
    }
  }

  // Perform a simulation dryrun step after setting simulation clock to given time
  doSimStepAt(currentTimeMs) {
    const currentDate = new Date(currentTimeMs);
    this.sim.setDate(currentDate); // Update sim clock, UI will read from it
    return this.sim.update();
  }

  // Narrows down a time bracket using the regula falsi method and returns a [x1, x2]
  // high-precision bracket. The order is preserved, if start < end, then t1 < t2.
  linearSolver(start, end, fStart, fEnd, tolerance, maxIterations = 10) {
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
        const fNext = this.doSimStepAt(xNext);

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
  findNextEvent(eventName, forward) {
    // Set the external function that will be called in Simulation
    // loop to pick the specific value according to event type
    const item = eventDefinitions.find(item => item.name === eventName);
    const fun = this.eventFunctions[item.type];
    this.sim.setDryRunFunction(fun);

    // Find a zero is not enough, we have to ensure the 'direction'
    // of the function, so to disambiguate between rise and set events
    // Transit value is always positive after transit (afternoon), and
    // negative before (morning), so it behaves like rise.
    const first_sign = (item.ascending === forward ? -1 : 1); // Sign of range's first start point

    // Chose a proper time step for findBracket.
    // Step should be large enough for efficency but also smaller than the minimum
    // interval between two consecutive (signed) zeros of the function.
    const period = (item.period * 24 * 3600 * 1000) * 0.75; // in msec
    const step = (forward ? period : -period);

    // Find initial [x1, x2] range
    const tolerance = 0.01;
    const start = this.sim.getTime() + step / 18; // Start some minutes before/after now
    const [x1, x2, f1, f2] = this.findBracket(start, step, first_sign, tolerance);

    // Find the zero of the function, that must be monotonic in [x1, x2], and return a
    // high-precision range [z1, z2] satisfying the condition (x1 < x2) == (z1 < z2)
    const [z1, z2] = this.linearSolver(x1, x2, f1, f2, tolerance);

    this.sim.setDryRunFunction(null);
    this.reset(); // Reset stale states
    return forward ? z2 : z1; // Pick the timestamp _after_ the event
  }
};
