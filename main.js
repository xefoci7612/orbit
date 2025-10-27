// ============================================================================
// UI CONTROLS
// ============================================================================

'use strict';

import Simulation, { VALIDATE } from './simulation.js';
import { CELESTIAL_EVENTS } from './events.js';
import { runValidation } from './validate.js';

const toDateStr = d => d.toISOString().slice(0, 10);
const toTimeStr = d => d.toISOString().slice(11, 11 + 5); // HH:mm

// Slider value moves linearly, we use a power law for actual speed
const SIM_TIME_SPEED_UP = 60; // 1 simulation hour elapses in 1 minute
const S0 = 30;
const K0 = 10 * 10000 / 9600;
const sliderToSpeed = s => K0 * Number(s) * Number(s) / (S0 * S0);
const speedToSlider = v => Math.sqrt(Math.abs(v / K0)) * S0;

// Formats a number like "+45.1°" or "-8.5°"
function formatDegrees(num) {
  const fixedNum = num.toFixed(1);
  const sign = num >= 0 ? '+' : '';
  return `${sign}${fixedNum}°`;
}

// Converts a decimal degree value into DMS (Degrees, Minutes, Seconds)
function toDMS(deg, isLatitude) {
  const tags = [['W', 'E'], ['S', 'N']];
  const dir = tags[isLatitude ? 1 : 0][deg >= 0 ? 1 : 0];
  const secs = Math.round(Math.abs(deg) * 3600);
  const [d, r] = [Math.floor(secs / 3600), (secs % 3600)];
  const [m, s] = [Math.floor(r / 60), (r % 60)];
  return `${d}° ${m}′ ${s}″ ${dir}`;
}

// Converts a DMS (Degrees, Minutes, Seconds) or a generic value string
// into a decimal degree value
function fromDMS(valueString) {
  let numericStr = valueString.trim();

  // Decimal only case (assume to be degrees)
  if (!numericStr.match(/[°'′"″NSEW]/i))
      return parseFloat(numericStr);

  // Determine sign and clean the string
  const sign = numericStr.match(/[SW]/i) ? -1 : 1;
  numericStr = numericStr.replace(/[NSEW]/i, '');

  // DMS parsing
  const matches = [];
  matches.push(numericStr.match(/([\d\.]+)\s*°/));    // Degrees
  matches.push(numericStr.match(/([\d\.]+)\s*['′]/)); // Minutes
  matches.push(numericStr.match(/([\d\.]+)\s*["″]/)); // Seconds

  // Accumulate the total seconds
  const k = [3600, 60, 1];
  const totalSeconds = matches.reduce((total, match, index) => {
    return total + (match ? parseFloat(match[1]) * k[index] : 0);
  }, 0);

  return (Math.abs(totalSeconds) / 3600) * sign; // secs to degrees
}

// Init simulation
const simulation = new Simulation(VALIDATE);
const renderer = simulation.getRenderer();

// Camera views
const slots = new Array(10).fill(null);

// DOM elements
const elementIds = {
  loaderOverlay: 'loader-overlay',
  crosshair: 'crosshair',
  simDate: 'sim-date',
  simTime: 'sim-time',
  btnReverse: 'btn-reverse',
  btnPlay: 'btn-play',
  btnReset: 'btn-reset',
  eventGroup: 'event-group',
  prevEvent: 'prev-event',
  eventIndicator: 'event-indicator',
  nextEvent: 'next-event',
  speedSlider: 'speed-slider',
  speedOut: 'speed-out',
  btnSpeed1x: 'btn-speed-1x',
  observerBtn: 'observer-btn',
  observerCam: 'observer-cam',
  satelliteBtn: 'satellite-btn',
  observerDataGroup: 'observer-data-group',
  azimuthOut: 'azimuth-out',
  elevationOut: 'elevation-out',
  latitudeOut: 'latitude-out',
  longitudeOut: 'longitude-out',
  satDataGroup: 'sat-data-group',
  satHeightOut: 'sat-height-out',
  satGrnSpeedOut: 'sat-grn-speed-out',
  satOrbitSpeedOut: 'sat-orbit-speed-out',
  lockIndicator: 'lock-indicator',
};

const elem = {};
Object.entries(elementIds).forEach(([key, id]) => {
  elem[key] = document.getElementById(id);
});

// A state-machine to handle observer UI transitions
class ObserverUI {
  constructor() {
    this.state = "idle";
  }

  // Method as arrow function automatically binds `this` to the ObserverUI class
  // object, this is needed for an event listener
  #placeObserverAt = (event) => {

    console.assert(this.state == "place", "ObserverUI: clicked while not placing!");

    // Pick the object clicked on by the user
    const hit = simulation.pickObject(event.clientX, event.clientY);
    if (!hit) {
      this.setState("idle");
      return;
    }
    // Create an observer view on the clicked point and associate it
    // to the observer camera
    const slotIdx = parseInt(elem.observerCam.dataset.slot);
    const viewIdx = simulation.createObserverView(hit.object, hit.point);
    slots[slotIdx] = viewIdx;
    this.setState("active");
  };

  reset() {
    if (this.state !== "idle") this.setState("idle");
  }

  setState(newState) {
    const curState = this.state;

    if (newState === "place") {

      console.assert(curState === "idle", "ObserverUI: place while not idle!");

      document.body.classList.add('observer-mode'); // change pointer
      elem.observerBtn.classList.add('enabled');
      renderer.domElement.addEventListener('click', this.#placeObserverAt);

    } else if (newState === "active") {

      console.assert(curState === "place", "ObserverUI: active while not placing!");

      document.body.classList.remove('observer-mode');
      elem.observerCam.classList.add('saved');
      elem.observerCam.hidden = false;
      renderer.domElement.removeEventListener('click', this.#placeObserverAt);

    } else if (newState === "idle") {

      console.assert(curState !== "idle", "ObserverUI: already idle!");

      document.body.classList.remove('observer-mode');
      elem.observerBtn.classList.remove('enabled');
      elem.observerCam.hidden = true;
      elem.eventGroup.hidden = true;
      renderer.domElement.removeEventListener('click', this.#placeObserverAt);
      if (curState === "active")
        simulation.disposeObserverView();
    }
    this.state = newState;
  }
};
const observerUI = new ObserverUI();

function exitSatelliteView() {
  simulation.exitSatelliteView();
  const satSlotIndex = parseInt(elem.satelliteBtn.dataset.slot);
  slots[satSlotIndex] = null;
}

// Simple object literal to handle long/short mouse clicks
const longPressHandler = {
  btnIdx: null,
  timer: null,

  onMouseDown(index) {
    this.btnIdx = index;
    this.timer = Date.now();
  },

  onMouseUp(index) {
    const sameBtn = this.btnIdx === index;
    this.btnIdx = null;
    const isLongClick = Date.now() - this.timer > 500;
    return sameBtn ? isLongClick : null;
  }
};

function camBtnPressed(slotIndex) {
  longPressHandler.onMouseDown(slotIndex);
}

function camBtnReleased(slotIndex) {
  const longClick = longPressHandler.onMouseUp(slotIndex);
  if (longClick === null)
    return; // released on a different button

  const isObserverBtn = this === elem.observerCam;
  const isSatelliteBtn = this === elem.satelliteBtn;

  // Long click: clear the slot
  if (longClick) {
    if (isObserverBtn || isSatelliteBtn) {
      return;
    }
    const viewIndex = slots[slotIndex];
    if (viewIndex === null) {
      return;
    }
    simulation.disposeView(viewIndex);
    slots[slotIndex] = null;
    this.classList.remove('saved');

  } else {
    // Normal click: switch the view
    let viewIndex = slots[slotIndex];

    // Satellite is a special case because we recreate/destroy
    // the view everytime we enter/exit
    if (isSatelliteBtn && !simulation.isSatelliteView()) {

      console.assert(viewIndex === null, "Satellite slot is not empty!");

      viewIndex = simulation.createSatelliteView();
      slots[slotIndex] = viewIndex;
    }

    // If slot is empty save the camera preset
    if (viewIndex === null) {

      // Saving a cam view while on observer or satellite is not supported
      if (simulation.isObserverView() || simulation.isSatelliteView()) {
        return;
      }

      console.assert(!isSatelliteBtn, "Satellite slot is empty!");
      console.assert(!isObserverBtn, "Observer slot is empty!");

      viewIndex = simulation.cloneView();
      slots[slotIndex] = viewIndex;
      this.classList.add('saved');
    }

    // We are leaving satellite view
    if (simulation.isSatelliteView() && !isSatelliteBtn) {
      exitSatelliteView();
    }

    // Finally switch to the new view
    simulation.setActiveView(viewIndex);
  }
  syncUI();
}

// Bind all camera buttons to click handlers
document.querySelectorAll('.cam-btn').forEach((b,i)=>{
  const slotIndex = parseInt(b.dataset.slot);
  b.onmousedown = camBtnPressed.bind(b, slotIndex);
  b.onmouseup = camBtnReleased.bind(b, slotIndex);
});

// Sync UI with simulation state
function syncUI() {

  // Date and time
  const currentTime = simulation.getTime();
  const date = new Date(currentTime);
  elem.simDate.value = toDateStr(date);
  elem.simTime.value = toTimeStr(date);

  // Speed slider
  const speed = simulation.speed();
  elem.speedSlider.value = speedToSlider(speed);
  if (speed !== 0) { // keep direction if speed == 0
    elem.btnReverse.textContent = speed < 0 ? '◅' : '▻';
  }
  const realSpeed = (SIM_TIME_SPEED_UP * Math.abs(speed)).toFixed(1);
  elem.speedOut.textContent = (speed < 0 ?'−':'') + realSpeed +' ×';

  // Lock icon
  if (simulation.isViewLocked()){
    elem.lockIndicator.classList.add('active');
  } else {
    elem.lockIndicator.classList.remove('active');
  }

  // Active camera button
  document.querySelectorAll('.cam-btn').forEach(btn => {
    btn.classList.remove('active');
  });
  const viewIdx = simulation.getActiveView();
  const slotIdx = slots.indexOf(viewIdx);
  const activeBtn = document.querySelector(`.cam-btn[data-slot="${slotIdx}"]`);
  if (activeBtn !== null)
    activeBtn.classList.add('active');

  // Crosshair in observer view
  if (simulation.isObserverView()) {
    crosshair.classList.add('visible');
  } else {
    crosshair.classList.remove('visible');
  }
}

function sceneDoubleClicked(mouseX, mouseY) {
  const hit = simulation.pickObject(mouseX, mouseY);
  if (!hit || simulation.isSatelliteView() || simulation.isObserverView())
    return;

  // Clicked again on the already locked object?
  if (simulation.isOrbitLocked(hit.object)) {
    simulation.unlockFromOrbit();
  } else {
    simulation.lockToOrbit(hit.object, hit.point);
  }
  syncUI();
};

function changeDateTime() {
  const newDate = new Date(elem.simDate.value + 'T' + elem.simTime.value + 'Z'); // UTC (Zulu time)
  simulation.setDate(newDate);
  syncUI();
}

elem.simDate.onchange = changeDateTime;
elem.simTime.onchange = changeDateTime;

elem.speedSlider.oninput = () => {
  const speed = sliderToSpeed(elem.speedSlider.value);
  const sign = elem.btnReverse.textContent === '◅' ? -1 : 1;
  simulation.setSpeed(sign * speed);
  syncUI();
};

elem.btnSpeed1x.onclick = () => {
  elem.speedSlider.value = speedToSlider(1 / SIM_TIME_SPEED_UP);
  elem.speedSlider.oninput();
}

elem.observerBtn.onclick = () => {
  const newState = observerUI.state === "idle" ? "place" : "idle";
  observerUI.setState(newState);
};

elem.btnReverse.onclick = () => {
  simulation.reverseSpeed();
  syncUI();
};

elem.btnPlay.onclick = () => {
  const isPlay = simulation.togglePause();
  elem.btnPlay.textContent = isPlay ? '⏸' : '▶';
  if (isPlay) {
    elem.eventGroup.hidden = true;
  }
};

elem.btnReset.onclick = () => {
  if (simulation.isSatelliteView())
    exitSatelliteView();
  observerUI.reset();
  simulation.reset();
  syncUI();
};

function updateEventDisplay(event) {
  elem.eventIndicator.classList.remove(...CELESTIAL_EVENTS.names);
  elem.eventIndicator.classList.add(event);
  elem.eventIndicator.textContent = elem.eventIndicator.dataset[event];
  elem.eventIndicator.title = event.charAt(0).toUpperCase() + event.slice(1);
}

simulation.on(CELESTIAL_EVENTS.event, (name) => {
  updateEventDisplay(name);
  elem.eventGroup.hidden = false;
  const inPause = (elem.btnPlay.textContent === '▶');
  if (!inPause)
    elem.btnPlay.onclick();
});

elem.prevEvent.onclick = () => {
  const eventName = elem.eventIndicator.classList[0];
  const timeForward = elem.btnReverse.textContent == '▻'; // speed can be 0
  simulation.goToNextEvent(eventName, false, timeForward);
  syncUI();
};

elem.nextEvent.onclick = () => {
  const eventName = elem.eventIndicator.classList[0];
  const timeForward = elem.btnReverse.textContent == '▻'; // speed can be 0
  simulation.goToNextEvent(eventName, true, timeForward);
  syncUI();
};

addEventListener('dblclick', (event) => {
  sceneDoubleClicked(event.clientX, event.clientY);
});

addEventListener('resize', () => {
  simulation.resize();
});

const setEditMode = (element, isFocused) => {
  element.setAttribute('data-manual-edit', isFocused.toString());
};

const onLatLonBlur = () => {
  // When the user clicks from lat to lon, the blur event fires, but the
  // focus on the new element hasn't happened yet, so use a small timeout
  setTimeout(() => {
    // Check if the currently focused element in the entire document is one of our inputs.
    const isFocusStillOnInputs =
      document.activeElement === elem.latitudeOut ||
      document.activeElement === elem.longitudeOut;

    // If focus has truly left the input group, reset edit mode
    // and commit the changes
    if (!isFocusStillOnInputs) {
      setEditMode(elem.latitudeOut, false);
      setEditMode(elem.longitudeOut, false);

      const latDeg = fromDMS(elem.latitudeOut.value);
      const lonDeg = fromDMS(elem.longitudeOut.value);

      if (!isNaN(latDeg) && !isNaN(lonDeg)) {
        simulation.placeObserverAt(latDeg, lonDeg);
      }
    }
  }, 0);
};

elem.latitudeOut.onfocus = () => setEditMode(elem.latitudeOut, true);
elem.longitudeOut.onfocus = () => setEditMode(elem.longitudeOut, true);
elem.latitudeOut.onblur = onLatLonBlur;
elem.longitudeOut.onblur = onLatLonBlur;

function animationLoop() {

  // Perform a simulation step and render
  const simData = simulation.update();

  elem.observerDataGroup.hidden = true;
  elem.satDataGroup.hidden = true;

  // If we are in observer/satellite view update dynamic info
  if (simData !== null) {
    if (simulation.isObserverView()) {
      elem.observerDataGroup.hidden = false;
      elem.azimuthOut.textContent = formatDegrees(simData.azimuth);
      elem.elevationOut.textContent = formatDegrees(simData.elevation);

      // Only update the UI if the user is NOT currently editing the fields
      const isEditingLat = elem.latitudeOut.getAttribute( 'data-manual-edit') === 'true';
      const isEditingLon = elem.longitudeOut.getAttribute('data-manual-edit') === 'true';
      if (!isEditingLat && !isEditingLon) {
        elem.latitudeOut.value  = toDMS(simData.latitude, true);
        elem.longitudeOut.value = toDMS(simData.longitude, false);
      }
    } else if (simulation.isSatelliteView()) {
      elem.satDataGroup.hidden = false;
      elem.satHeightOut.textContent   = `${simData.height} Km`;
      elem.satGrnSpeedOut.textContent = `${simData.groundSpeed} Km/h`;
      elem.satOrbitSpeedOut.textContent = `${simData.inertialSpeed} Km/h`;
    }
  }

  // Schedule next step
  requestAnimationFrame(animationLoop);
}

window.addEventListener('load', () => {
  // Remove loder overlay
  elem.loaderOverlay.classList.add('hidden'); // trigger fade-out transition
  elem.loaderOverlay.addEventListener('transitionend', () => {
    elem.loaderOverlay.remove();
  });

  // Sync UI with simulation initial state
  simulation.update();
  elem.btnReset.onclick();

  // If in validation mode run sampling code and exit
  if (VALIDATE) {
    runValidation(simulation);
    return;
  }

  // Start animation loop after all is loaded
  animationLoop();
});
