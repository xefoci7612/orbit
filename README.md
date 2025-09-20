# Interactive Earth-Moon System Simulation

This is a real-time 3D visualization of the Earth-Moon system, allowing you to explore the celestial mechanics of our planet and its natural satellite. You can control the flow of time, observe the system from any angle, and even place a virtual observer on the surface of the Earth or Moon.

## Main Features

*   **Time Manipulation**: Control the pace and direction of the simulation. You can play, pause, reverse time, and adjust the simulation speed from a slow crawl to many times faster than real-time.
*   **Accurate Orbital Mechanics**: The simulation uses precise orbital parameters from JPL Horizons for September 10, 2025, to accurately model the Moon's elliptical and inclined orbit, as well as its precession.
*   **Interactive 3D Camera**: Freely navigate the scene by rotating (left-click and drag), panning (right-click and drag), and zooming (mouse wheel).
*   **Camera Presets**: Save and recall up to three custom camera positions. A fourth slot is reserved for the special "Observer" view.
*   **Lock-on View**: Lock the camera onto either the Earth or the Moon. The camera will maintain its position relative to the moving and rotating body, allowing you to follow it through its cosmic journey.
*   **Observer Mode**: Place an observer anywhere on the surface of the Earth or Moon to see the sky from that perspective. When an observer is placed, the simulation will automatically pause at key astronomical events like moonrise and moonset.

## How to Use the Simulation

### The Control Panel (HUD)

The controls are located in a semi-transparent panel at the bottom of the screen.

#### **Time Controls**

*   **Date and Time Input**: On the left, you can set a specific date and time to jump to that moment in the simulation.
*   **Elapsed Time**: Next to the date, a counter shows how many days have passed since the simulation's start date of September 10, 2025.
*   **Playback Buttons (Center)**:
    *   `▻` / `◅`: Toggles the direction of time (forward or reverse).
    *   `⏸` / `▶`: Pauses or resumes the simulation.
    *   `⟲`: Resets the simulation to its initial date, time, and camera view.
*   **Speed Slider (Right)**:
    *   Drag the slider to increase or decrease the simulation speed. The value shown is a multiplier of the base speed (where 1× represents one simulation hour passing in one real-time minute).

#### **Camera Controls & Modes**

A second control row appears for managing camera views.

*   **Saving a View**: Navigate to a desired viewpoint and then click an empty camera button (`1`, `2`, or `3`) to save it. The button will light up, indicating a view is saved to that slot.
*   **Recalling a View**: Click a saved camera button to instantly return to that viewpoint.
*   **Clearing a View**: Long-press (click and hold for about half a second) a saved camera button to clear its stored position.

#### **Special Views**

*   **Lock-on Feature**:
    1.  **Double-click** on either the Earth or the Moon.
    2.  The camera will lock onto the selected body, automatically following its rotation and orbit. A lock icon (`🔒`) will appear in the camera panel to indicate it is active.
    3.  A green marker will appear on the surface where you clicked.
    4.  To unlock, simply double-click the same celestial body again.

*   **Observer Mode**:
    1.  Click the pin icon (`📍`) to enter Observer Mode. The mouse cursor will change to a crosshair.
    2.  Click on the surface of the Earth or Moon where you'd like to place your observer.
    3.  The view will instantly change to the perspective of that observer, looking out from the surface.
    4.  The button `4` becomes active as the dedicated slot for this view.
    5.  The simulation will now automatically pause when the Moon rises or sets from the observer's location.
    6.  To exit Observer Mode, click the pin icon (`📍`) again.
