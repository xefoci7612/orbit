// ============================================================================
// UI CONTROLS
// ============================================================================

'use strict';

import Simulation from './simulation.js';

const toDateStr = d => d.toISOString().slice(0, 10);
const toTimeStr = d => d.toISOString().slice(11, 11 + 5); // HH:mm

// Slider value moves linearly, we use a power law for actual speed
const sliderToSpeed = v => Number(v) * Number(v) / 2.5;
const speedToSlider = v => Math.sqrt(Math.abs(v * 2.5));

const simulation = new Simulation();
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

    const hit = simulation.getClickedCandidate(event.clientX, event.clientY);
    if (!hit)
      return;

    // Place observer camera
    const slotIdx = parseInt(elem.observerCam.dataset.slot) - 1;
    const viewIdx = simulation.setObserver(hit.object, hit.point);
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
      simulation.removeObserver();

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
  const [date, elapsedMsec] = simulation.getDate();
  elem.simDate.value = toDateStr(date);
  elem.simTime.value = toTimeStr(date);
  const elapsedDays = elapsedMsec / (24 * 3600 * 1000);
  elem.simElapsed.textContent = (elapsedDays).toFixed(1)+' d';
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
  const hit = simulation.getClickedCandidate(mouseX, mouseY);
  if (!hit)
    return;
  // Clicked again on the already locked object?
  if (simulation.isLocked(hit.object)) {
    simulation.unlockCamera();
    elem.linkIndicator.classList.remove('active');
  } else {
    simulation.lockCameraTo(hit.object, hit.point);
    elem.linkIndicator.classList.add('active');
  }
};

elem.simDate.onchange = dateTimeChanged;
elem.simTime.onchange = dateTimeChanged;
elem.speedSlider.oninput = speedSliderChanged;

elem.observerBtn.onclick = () => {
  const newState = observerUI.state === "idle" ? "place" : "idle";
  observerUI.setState(newState);
};

elem.btnPlay.onclick = ()=>{
  elem.btnPlay.textContent = simulation.togglePause() ? '⏸' : '▶';
};

elem.btnReverse.onclick = ()=> {
  simulation.reverseSpeed();
  updateSpeedUI();
};

elem.btnReset.onclick = ()=>{
  simulation.reset();
  updateSpeedUI();
  updateDateTimeUI();
};

addEventListener('dblclick', (event) => {
  sceneDoubleClicked(event.clientX, event.clientY);
});

addEventListener('resize', () => {
  simulation.resize();
});

function animationLoop() {

  // Perform a simulation step and render
  const events = simulation.update();

  // Events detected during simulation step
  if (events.atRise || events.atSet)
      elem.btnPlay.onclick();

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

  // Start animation loop after all is loaded
  simulation.setDate(new Date()); // now
  animationLoop();
});
