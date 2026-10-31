# Integration checks for agents

Checks to run after porting TimelineFX to a renderer or engine, or after changing the integration. You can't watch the
window, so each one prints something to compare. The porting checklist in implementations/README.md says what to build;
this says how to tell it works.

**First:** route the library's messages to a log you can read with `tfx_SetPrintCallback`, if the host has no console,
and read every line starting `TimelineFX tip:`. The tip saying an effect "spawns at user spawn locations" is expected
once per template when you add locations straight after. Act on every other one.

## 1. The stage is healthy

Run the scene for a few seconds at the update rate the game will use, then print `tfx_DescribeStage` straight after a
frame's `tfx_GetGPUParticleFrame` (see gpu-particles.md for the call). Done when it has only its counts and GPU lines: no `Blocked:`, `Dropping spawns:`,
`Unfetched:` or `Library ...` lines. A `No spawn locations right now` line is fine only if that effect is meant to be
idle at that moment.

## 2. The GPU matches the CPU

Run the same scene twice, with the same `tfx_SetStageSeed` and the same schedule of effects and locations: once as is,
and once with `tfx_SetStageGPUSimulation(stage, false)` set before the effects are added. Average over a few seconds:

- CPU run: `tfx_GetParticleCount`.
- GPU run: `tfx_GetParticleCount` (the emitters still on the CPU) + `tfx_GetGPUParticleStats(stage).estimated_live_particles`.

Done when the two are within 10% (they're usually within 3%). Never compare against `slots_in_use`, which counts dead
particles still holding their slots. A big gap with no dropped spawns or blocked emitters means the GPU integration is
losing particles: check that every fetch is dispatched exactly once, and that each layer draws its GPU range.

## 3. No drops at the real update rate

`tfx_GetGPUParticleStats(stage).dropped_spawns` stays 0 and `blocked_emitters` stays 0 at the game's update rate, with
the most locations the game will ever have playing at once. Ring sizes assume about 60 updates a second, so test at the
real rate. Size with `tfx_SetTemplateGPUCapacityForLocations` (gpu-particles.md).

## 4. Fetches keep up

`tfx_GetStageDiagnostics(stage).unfetched_gpu_updates` is 0 straight after each frame's fetch, and never builds up
across frames. If it climbs, the host is updating the stage on frames it doesn't fetch and dispatch.

## 5. The update rate doesn't show

Run at 20 and at 400 updates a second. Particles should move and look the same as at 60 (the vertex shader
interpolates between updates). If they stutter or jump, the instance buffers aren't swapping exactly once per update,
or a host that double buffers its own instances is breaking a rule in "Engines with their own instance buffers" in
implementations/README.md. At 20 updates a second, GPU rings sized for 60 drop spawns, so check 3 only applies at the
real rate.
