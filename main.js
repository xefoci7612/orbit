// ============================================================================
// UI CONTROLS
// ============================================================================

import Simulation from './simulation.js';

const simulation = new Simulation();
const renderer = simulation.getRenderer();

const dateInp  = document.getElementById('simDate');
const timeInp  = document.getElementById('simTime');
const elapsedS = document.getElementById('simElapsed');
const playBtn  = document.getElementById('btnPlay');
const revBtn   = document.getElementById('btnReverse');
const resBtn   = document.getElementById('btnReset');
const speedSl  = document.getElementById('speedSlider');
const speedOut = document.getElementById('speedOut');
const camRow   = document.getElementById('camRow');
const camStore = document.getElementById('camStore');
const camClose = document.getElementById('camClose');
const linkIcon = document.getElementById('linkIndicator');
const observerBtn = document.getElementById('observerBtn');
const observerCam = document.getElementById('observerCam');

const toDateStr = d => d.toISOString().slice(0, 10);
const toTimeStr = d => d.toISOString().slice(11, 19);

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

// Slider value moves linearly, we use a power law for actual speed
const sliderToSpeed = v => Number(v) * Number(v) / 2.5;
const speedToSlider = v => Math.sqrt(Math.abs(v * 2.5));

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
  if (Lock.active && hit.object === Lock.lockedObj) {
    simulation.unlockCamera();
    linkIcon.classList.remove('active');
  } else {
    simulation.lockCameraTo(hit.object, hit.point);
    linkIcon.classList.add('active');
  }
};

function observerTargetClicked(mouseX, mouseY) {
  if (!document.body.classList.contains('observer-mode'))
    return;

  const hit = simulation.getClickedCandidate(mouseX, mouseY);
  if (!hit)
    return;

  // Setup observer camera: use the 1-based data-slot value, converted to a 0-based index
  const slotIdx = parseInt(observerCam.dataset.slot) - 1;
  const viewIdx = simulation.setObserver(hit.object, hit.point);
  slots[slotIdx] = viewIdx;

  document.body.classList.remove('observer-mode');
  renderer.domElement.removeEventListener('click', observerClickListener);
  observerClickListener = null;
  observerCam.classList.add('saved');
  observerCam.hidden = false;
}

dateInp.onchange = dateTimeChanged;
timeInp.onchange = dateTimeChanged;
speedSl.oninput = speedSliderChanged;

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

let observerClickListener = null;

observerBtn.onclick = () => {
  const active = simulation.observerActive() || observerClickListener;
  if (!active) {
    document.body.classList.add('observer-mode');
    observerBtn.style.background = 'rgba(255,255,255,.5)'; // Highlight the button

    // Add temporary click listener for placing observer, cannot be permanent
    // otherwise would raise even when pressing this button
    observerClickListener = (event) => { observerTargetClicked(event.clientX, event.clientY); };
    renderer.domElement.addEventListener('click', observerClickListener);
  } else { // Exit observer mode
    document.body.classList.remove('observer-mode');
    observerBtn.style.background = ''; // Reset button style
    observerCam.hidden = true;

    if (observerClickListener) {
      renderer.domElement.removeEventListener('click', observerClickListener);
      observerClickListener = null;
    }
    simulation.removeObserver();
  }
};

const slots = new Array(4).fill(null); // Will store view indices from ViewManager
let pressedSlotIndex = null;
let clickTimer = null;

function CamBtnPressed(slotIndex) {
  pressedSlotIndex = slotIndex;
  clickTimer = Date.now(); // start timer for double-click detection
}

function CamBtnReleased(slotIndex) {
  const isObserverCamera = (this === observerCam);
  const longClick = (Date.now() - clickTimer) > 400; // 400 ms threshold
  const sameBtn = pressedSlotIndex === slotIndex;
  pressedSlotIndex = null; // reset
  clickTimer = null;

  if (!sameBtn)
    return;

  // Long click: clear the slot
  if (longClick && !isObserverCamera) {
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

// Bind all camera buttons to click and double-click handlers
document.querySelectorAll('.camBtn').forEach((b,i)=>{
  // Use the data-slot value (1-based) converted to a 0-based index
  const slotIndex = parseInt(b.dataset.slot) - 1;
  b.onmousedown = CamBtnPressed.bind(b, slotIndex);
  b.onmouseup = CamBtnReleased.bind(b, slotIndex);
});

addEventListener('dblclick', (event) => {
  if (!document.body.classList.contains('observer-mode'))
    sceneDoubleClicked(event.clientX, event.clientY);
});

addEventListener('resize', () => {
  simulation.resize();
});

function animationLoop() {
  simulation.update();
  updateDateTimeUI();
  requestAnimationFrame(animationLoop);
}

// Finally start animation
animationLoop();
