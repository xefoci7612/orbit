// ============================================================================
// UI CONTROLS
// ============================================================================

import Simulation from './simulation.js';

const simulation = new Simulation();
const renderer = simulation.getRenderer();

// View indices from ViewManager
const slots = new Array(10).fill(null);

const dateInp  = document.getElementById('simDate');
const timeInp  = document.getElementById('simTime');
const elapsedS = document.getElementById('simElapsed');
const playBtn  = document.getElementById('btnPlay');
const revBtn   = document.getElementById('btnReverse');
const resBtn   = document.getElementById('btnReset');
const speedSl  = document.getElementById('speedSlider');
const speedOut = document.getElementById('speedOut');
const linkIcon = document.getElementById('linkIndicator');
const observerBtn = document.getElementById('observerBtn');
const observerCam = document.getElementById('observerCam');
const loaderOvr = document.getElementById('loader-overlay');

const toDateStr = d => d.toISOString().slice(0, 10);
const toTimeStr = d => d.toISOString().slice(11, 11 + 5); // HH:mm

// Slider value moves linearly, we use a power law for actual speed
const sliderToSpeed = v => Number(v) * Number(v) / 2.5;
const speedToSlider = v => Math.sqrt(Math.abs(v * 2.5));

// Simple class to handle long/short mouse clicks
class LongPressHandler {
  constructor() {
    this.slotIdx = null;
    this.clickTimer = null;
  }
  onMouseDown(index) {
    this.slotIdx = index;
    this.clickTimer = Date.now(); // start timer for double-click detection
  }
  onMouseUp(index) {
    const sameBtn = (this.slotIdx == index);
    this.slotIdx = null;
    const isLongClick = (Date.now() - this.clickTimer) > 500; // 500 ms
    return !sameBtn ? null : isLongClick;
  }
};
const longPressHandler = new LongPressHandler;

// A state-machine to handle observer UI transitions
class ObserverUI {
  constructor() {
    this.state = "idle";
    this.onClick = null;
  }
  getState() {
    return this.state;
  }
  setState(state) {
    if (state == "place") {
      console.assert(this.state == "idle", "ObserverUI: place while not idle!");
      document.body.classList.add('observer-mode');
      observerBtn.style.background = 'rgba(255,255,255,.5)'; // Highlight the button

      // Add temporary click listener for placing observer, cannot be permanent
      // otherwise would raise even when pressing this button
      this.onClick = (event) => { observerTargetClicked(event.clientX, event.clientY); };
      renderer.domElement.addEventListener('click', this.onClick);
    }
    else if (state == "active") {
      console.assert(this.state == "place", "ObserverUI: active while not placing!");
      document.body.classList.remove('observer-mode');
      observerCam.classList.add('saved');
      observerCam.hidden = false;

      renderer.domElement.removeEventListener('click', this.onClick);
      this.onClick = null;
    }
    else if (state == "idle") {
      console.assert(this.state !== "idle", "ObserverUI: already idle!");
      document.body.classList.remove('observer-mode');
      observerBtn.style.background = '';
      observerCam.hidden = true;

      // User clicked observer button while still looking for a place
      if (state == "place") {
        renderer.domElement.removeEventListener('click', this.onClick);
        this.onClick = null;
      }
    }
    this.state = state;
  }
};
const observerUI = new ObserverUI();

function updateDateTimeUI() {
  const [date, elapsedMsec] = simulation.getDate();
  dateInp.value = toDateStr(date);
  timeInp.value = toTimeStr(date);
  const elapsedDays = elapsedMsec / (24 * 3600 * 1000);
  elapsedS.textContent = (elapsedDays).toFixed(1)+' d';
}

function dateTimeChanged() {
  const newDate = new Date(dateInp.value + 'T' + timeInp.value + 'Z');
  simulation.setDate(newDate);
  updateDateTimeUI();
}

function updateSpeedUI() {
  const speed = simulation.speed();
  speedSl.value = speedToSlider(speed);
  if (speed !== 0) // keep direction if speed == 0
    revBtn.textContent = speed < 0 ? '◅' : '▻';
  speedOut.textContent = (speed < 0 ?'−':'') + Math.abs(speed).toFixed(1)+' ×';
}

function speedSliderChanged() {
  const speed = sliderToSpeed(speedSl.value);
  const sign = revBtn.textContent === '◅' ? -1 : 1;
  simulation.setSpeed(sign * speed);
  updateSpeedUI();
};

function sceneDoubleClicked(mouseX, mouseY) {
  const hit = simulation.getClickedCandidate(mouseX, mouseY);
  if (!hit)
    return;
  // Clicked again on the already locked object?
  if (simulation.isLocked(hit.object)) {
    simulation.unlockCamera();
    linkIcon.classList.remove('active');
  } else {
    simulation.lockCameraTo(hit.object, hit.point);
    linkIcon.classList.add('active');
  }
};

function observerTargetClicked(mouseX, mouseY) {
  console.assert(observerUI.getState() == "place", "ObserverUI: clicked while not placing!");
  const hit = simulation.getClickedCandidate(mouseX, mouseY);
  if (!hit)
    return;
  // Setup observer camera: use the 1-based data-slot value, converted to a 0-based index
  const slotIdx = parseInt(observerCam.dataset.slot) - 1;
  const viewIdx = simulation.setObserver(hit.object, hit.point);
  slots[slotIdx] = viewIdx;
  observerUI.setState("active"); // update UI elements
}

playBtn.onclick = ()=>{
  playBtn.textContent = simulation.togglePause() ? '⏸' : '▶';
};

revBtn.onclick = ()=> {
  simulation.reverseSpeed();
  updateSpeedUI();
};

resBtn.onclick = ()=>{
  simulation.reset();
  updateDateTimeUI();
  updateSpeedUI();
};

observerBtn.onclick = () => {
  if (observerUI.getState() == "idle") {
    observerUI.setState("place");
  } else {
    observerUI.setState("idle");
    simulation.removeObserver();
  }
};

function CamBtnPressed(slotIndex) {
  longPressHandler.onMouseDown(slotIndex);
}

function CamBtnReleased(slotIndex) {
  const longClick = longPressHandler.onMouseUp(slotIndex);
  if (longClick === null)
    return;
  // Long click: clear the slot
  if (longClick && this !== observerCam) {
    slots[slotIndex] = null;
    this.classList.remove('saved');
  } else {
    // Normal click: either save or restore camera
    // according if slot is empty or not
    const viewIndex = slots[slotIndex];
    if (viewIndex !== null) {
      simulation.setActiveView(viewIndex);
    } else {
       // Save new view: create a new camera and controls
      const newViewIndex = simulation.cloneView();
      slots[slotIndex] = newViewIndex;
      this.classList.add('saved');
    }
  }
}

dateInp.onchange = dateTimeChanged;
timeInp.onchange = dateTimeChanged;
speedSl.oninput = speedSliderChanged;

// Bind all camera buttons to click and double-click handlers
document.querySelectorAll('.camBtn').forEach((b,i)=>{
  // Use the data-slot value (1-based) converted to a 0-based index
  const slotIndex = parseInt(b.dataset.slot) - 1;
  b.onmousedown = CamBtnPressed.bind(b, slotIndex);
  b.onmouseup = CamBtnReleased.bind(b, slotIndex);
});

addEventListener('dblclick', (event) => {
  if (observerUI.getState() != "place")
    sceneDoubleClicked(event.clientX, event.clientY);
});

addEventListener('resize', () => {
  simulation.resize();
});

function animationLoop() {
  const events = simulation.update();
  if (events.atRise || events.atSet)
      playBtn.onclick();
  updateDateTimeUI();
  requestAnimationFrame(animationLoop);
}

window.addEventListener('load', () => {
  // Remove loder overlay
  loaderOvr.classList.add('hidden'); // trigger fade-out transition
  loaderOvr.addEventListener('transitionend', () => {
    loaderOvr.remove();
  });

  // Start animation loop
  simulation.setDate(new Date()); // now
  animationLoop();
});
