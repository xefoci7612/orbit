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
const sliderToSpeed = v => Number(v) * Number(v) / 2.5;
const speedToSlider = v => Math.sqrt(Math.abs(v * 2.5));

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

const simulation = new Simulation(VALIDATE);
const renderer = simulation.getRenderer();

// View indices from ViewManager
const slots = new Array(10).fill(null);

// DOM elements
const elementIds = {
  simDate: 'sim-date',
  simTime: 'sim-time',
  simElapsed: 'sim-elapsed',
  btnPlay: 'btn-play',
  btnReverse: 'btn-reverse',
  btnReset: 'btn-reset',
  speedSlider: 'speed-slider',
  speedOut: 'speed-out',
  linkIndicator: 'link-indicator',
  observerBtn: 'observer-btn',
  observerCam: 'observer-cam',
  loaderOverlay: 'loader-overlay',
  observerDataGroup: 'observer-data-group',
  azimuthOut: 'azimuth-out',
  elevationOut: 'elevation-out',
  latitudeOut: 'latitude-out',
  longitudeOut: 'longitude-out',
  eventGroup: 'event-group',
  prevEvent: 'prev-event',
  eventIndicator: 'event-indicator',
  nextEvent: 'next-event',
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

  // Method as arrow function automatically binds `this` to the class object
  // Needed, we will use this function as event listener
  #placeObserverAt = (event) => {

    console.assert(this.state == "place", "ObserverUI: clicked while not placing!");

    const hit = simulation.pickObject(event.clientX, event.clientY);
    if (!hit)
      return;

    // Place observer camera, point is in world coordinates
    const slotIdx = parseInt(elem.observerCam.dataset.slot) - 1;
    const viewIdx = simulation.enterObserverView(hit.object, hit.point);
    slots[slotIdx] = viewIdx;
    this.setState("active");
  };

  setState(newState) {
    const curState = this.state;

    if (newState === "place") {

      console.assert(curState === "idle", "ObserverUI: place while not idle!");

      document.body.classList.add('observer-mode');
      elem.observerBtn.style.background = 'rgba(255,255,255,.5)';
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
      elem.observerBtn.style.background = '';
      elem.observerCam.hidden = true;
      simulation.exitObserverView();

      // User clicked observer button while still looking for a place
      if (curState === "place") {
        renderer.domElement.removeEventListener('click', this.#placeObserverAt);
      }
    }
    this.state = newState;
  }
};
const observerUI = new ObserverUI();

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

  // Long click: clear the slot
  if (longClick && this !== elem.observerCam) {
    const viewIndex = slots[slotIndex];
    if (viewIndex !== null) {
      simulation.disposeView(viewIndex);
      slots[slotIndex] = null;
      this.classList.remove('saved');
    }
  } else {
    // Normal click: either save or restore camera
    // according if slot is empty or not
    const viewIndex = slots[slotIndex];
    if (viewIndex === null) {
      slots[slotIndex] = simulation.cloneView();
      this.classList.add('saved');
    } else {
      simulation.setActiveView(viewIndex);
    }
  }
}

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

// Bind all camera buttons to click handlers
document.querySelectorAll('.cam-btn').forEach((b,i)=>{
  const slotIndex = parseInt(b.dataset.slot) - 1;
  b.onmousedown = camBtnPressed.bind(b, slotIndex);
  b.onmouseup = camBtnReleased.bind(b, slotIndex);
});

function dateTimeChanged() {
  const newDate = new Date(elem.simDate.value + 'T' + elem.simTime.value + 'Z'); // UTC (Zulu time)
  simulation.setDate(newDate);
  updateDateTimeUI();
}

function updateDateTimeUI() {
  const currentTime = simulation.getTime();
  const date = new Date(currentTime);
  elem.simDate.value = toDateStr(date);
  elem.simTime.value = toTimeStr(date);
  //const elapsedDays = elapsedMsec / (24 * 3600 * 1000);
  //elem.simElapsed.textContent = (elapsedDays).toFixed(1)+' d';
}

function speedSliderChanged() {
  const speed = sliderToSpeed(elem.speedSlider.value);
  const sign = elem.btnReverse.textContent === '◅' ? -1 : 1;
  simulation.setSpeed(sign * speed);
  updateSpeedUI();
};

function updateSpeedUI() {
  const speed = simulation.speed();
  elem.speedSlider.value = speedToSlider(speed);
  if (speed !== 0) // keep direction if speed == 0
    elem.btnReverse.textContent = speed < 0 ? '◅' : '▻';
  elem.speedOut.textContent = (speed < 0 ?'−':'') + Math.abs(speed).toFixed(1)+' ×';
}

function sceneDoubleClicked(mouseX, mouseY) {
  const hit = simulation.pickObject(mouseX, mouseY);
  if (!hit)
    return;

  // Clicked again on the already locked object?
  if (simulation.isOrbitLocked(hit.object)) {
    simulation.unlockCamera();
    elem.linkIndicator.classList.remove('active');
  } else if (!simulation.isObserverView()) {
    simulation.lockToOrbit(hit.object, hit.point);
    elem.linkIndicator.classList.add('active');
  }
};

elem.simDate.onchange = dateTimeChanged;
elem.simTime.onchange = dateTimeChanged;
elem.speedSlider.oninput = speedSliderChanged;

elem.observerBtn.onclick = () => {
  const newState = observerUI.state === "idle" ? "place" : "idle";
  observerUI.setState(newState);
  if (newState === "idle") elem.eventGroup.hidden = true;
};

elem.btnPlay.onclick = () => {
  const isPlay = simulation.togglePause();
  elem.btnPlay.textContent = isPlay ? '⏸' : '▶';
  if (isPlay) elem.eventGroup.hidden = true;
};

elem.btnReverse.onclick = () => {
  simulation.reverseSpeed();
  updateSpeedUI();
};

elem.btnReset.onclick = () => {
  simulation.reset();
  updateSpeedUI();
  updateDateTimeUI();
};

function updateEventDisplay(event) {
  elem.eventIndicator.classList.remove(...CELESTIAL_EVENTS.names);
  elem.eventIndicator.classList.add(event);
  elem.eventIndicator.textContent = elem.eventIndicator.dataset[event];
  elem.eventIndicator.title = event.charAt(0).toUpperCase() + event.slice(1);
}

function getNextEvent(delta) {
  const num = CELESTIAL_EVENTS.names.length;
  const curEvent = elem.eventIndicator.classList[0];
  const curId = CELESTIAL_EVENTS.names.indexOf(curEvent);
  const nextIndex = (curId + delta + num) % num; // in [0, 3] range
  return CELESTIAL_EVENTS.names[nextIndex];
}

elem.prevEvent.onclick = () => {
  const eventName = elem.eventIndicator.classList[0];
  simulation.findNextEvent(eventName, false);
  updateDateTimeUI();
};

elem.nextEvent.onclick = () => {
  const eventName = elem.eventIndicator.classList[0];
  simulation.findNextEvent(eventName, true);
  updateDateTimeUI();
};

addEventListener('dblclick', (event) => {
  sceneDoubleClicked(event.clientX, event.clientY);
});

addEventListener('resize', () => {
  simulation.resize();
});

simulation.on(CELESTIAL_EVENTS.event, (name) => {
  updateEventDisplay(name);
  elem.eventGroup.hidden = false;
  const inPause = (elem.btnPlay.textContent === '▶');
  if (!inPause)
    elem.btnPlay.onclick();
});

function animationLoop() {

  // Perform a simulation step and render
  const events = simulation.update();

  // If we are in observer view show/set coordinates
  if (events.azEl.length > 0) {
    elem.observerDataGroup.hidden = false;
    elem.azimuthOut.textContent = formatDegrees(events.azEl[0]);
    elem.elevationOut.textContent = formatDegrees(events.azEl[1]);

    // Only update the UI if the user is NOT currently editing the fields
    const isEditingLat = elem.latitudeOut.getAttribute( 'data-manual-edit') === 'true';
    const isEditingLon = elem.longitudeOut.getAttribute('data-manual-edit') === 'true';
    if (!isEditingLat && !isEditingLon) {
      elem.latitudeOut.value  = toDMS(events.latLon[0], true);
      elem.longitudeOut.value = toDMS(events.latLon[1], false);
    }
  } else {
    elem.observerDataGroup.hidden = true;
  }

  // Sync UI elements to new state
  updateDateTimeUI();

  // Schedule next step
  requestAnimationFrame(animationLoop);
}

window.addEventListener('load', () => {
  // Remove loder overlay
  elem.loaderOverlay.classList.add('hidden'); // trigger fade-out transition
  elem.loaderOverlay.addEventListener('transitionend', () => {
    elem.loaderOverlay.remove();
  });

  // If in validation mode run sampling code and exit
  if (VALIDATE) {
    runValidation(simulation);
    return;
  }

  // Start animation loop after all is loaded
  animationLoop();
});
