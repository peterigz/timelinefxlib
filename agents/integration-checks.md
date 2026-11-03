# Integration checks for agents

Checks to run after porting TimelineFX to a renderer or engine, or after changing the integration. You can't watch the
window, so each one prints something to compare. The porting checklist in implementations/README.md says what to build;
this says how to tell it works.

**First:** route the library's messages to a log you can read with `tfx_SetPrintCallback`, if the host has no console,
and read every line starting `TimelineFX tip:`, `TimelineFX error:` or `TimelineFX assert failed:`. The tip saying an effect "spawns at user spawn locations" is expected
once per effect when you add locations straight after. Act on every other one.

## 1. The stage is healthy

Run the scene for a few seconds at the update rate the game will use, then print `tfx_DescribeStage` straight after a
frame's `tfx_GetGPUParticleFrame` (see gpu-particles.md for the call). Done when it has only its counts and GPU lines: no `Blocked:`, `Dropping spawns:`,
`Unfetched:` or `Library ...` lines. A `No spawn locations right now` line is fine only if that effect is meant to be
idle at that moment. A `No spawn locations left` line means effects that have finished their bursts are piling up:
expire them, or add their locations with `tfxSpawnLocationAdd_expire_effect`. A library saved by an older editor loads and plays fine but shows a `Library ...` line with
`file_version_out_of_date` until it's re-saved in the editor.

## 2. The GPU matches the CPU

Run the same scene twice, with the same `tfx_SetStageSeed` and the same schedule of effects and locations: once as is,
and once with `tfx_SetStageGPUSimulation(stage, false)` set before the effects are added. Average over a few seconds:

- CPU run: `tfx_GetParticleCount`.
- GPU run: `tfx_GetParticleCount` (the emitters still on the CPU) + `tfx_GetGPUParticleStats(stage).estimated_live_particles`.

Done when the two are within 10%. The CPU count reads high on emitters whose lives vary: `tfx_GetParticleCount` includes
particles that died but still hold their place until the oldest in their emitter dies, and the estimate counts only the
living, so a gap of about half the life variation over the average life is expected (6% for lives of 2000 to 2250 ms).
Never compare against `slots_in_use`, which counts dead particles still holding their slots. To find which effect a gap
comes from, play one at a time, or read `tfx_GetEffectGPUParticleStats(stage, effect).estimated_live_particles` for
each. A big gap with no dropped spawns or blocked emitters means the GPU integration is losing particles: check that
every fetch is dispatched exactly once, and that every GPU range is drawn: one per layer, or one per key from
`tfx_GetNextInstanceKeyRange` when you use sort keys.

## 3. No drops at the real update rate

`tfx_GetGPUParticleStats(stage).dropped_spawns` stays 0 and `blocked_emitters` stays 0 at the game's update rate, with
the most locations the game will ever have playing at once, and through the catch-up frames a slow frame causes. Ring
sizes assume 60 updates a second unless `tfx_SetGPUSizingUpdateRate` was set, so test at the real rate. Size with
`tfx_SetTemplateGPUCapacityForLocations` and `tfx_SetGPUSizingUpdateRate` (gpu-particles.md). The count is since the stage was
created or last cleared, so `tfx_ClearStage` between test runs. `tfx_GetEffectGPUParticleStats` shows which effect drops.

## 4. Fetches keep up

`tfx_GetStageDiagnostics(stage).unfetched_gpu_updates` is 0 straight after each frame's fetch, and never builds up
across frames. If it climbs, the host is updating the stage on frames it doesn't fetch and dispatch.

## 5. The update rate doesn't show

Run at 20 and at 400 updates a second. Particles should move and look the same as at 60 (the vertex shader
interpolates between updates). If they stutter or jump, the instance buffers aren't swapping exactly once per update,
or a host that double buffers its own instances is breaking a rule in "Engines with their own instance buffers" in
implementations/README.md. At 20 updates a second, the rings of continuous GPU emitters sized for 60 can drop spawns:
each update spawns three times as many, and an update's slots are held until its longest lived particle dies. Single
shot bursts spawn in one update at any rate. So check 3 only applies at the real rate, unless
`tfx_SetGPUSizingUpdateRate(20)` was set first, in which case 20 must not drop either.

## 6. Libraries share one set of buffers

With more than one library on a stage that has `global_library_indexes` on (see "Drawing effects from more than one
library" in timelinefx.h):

- `tfx_UpdateGlobalGPUData` returns true after every load, and `tfx_DescribeStage` shows `Library indexes: global`
  with no line saying libraries didn't fit.
- Every effect added from each library gets a valid id (`tfx_EffectIDIsValid`). An invalid one means it was added before
  a sync placed its library.
- Print the four versions the backend uploads against (`tfx_GetGlobalParticlePropertiesVersion`,
  `tfx_GetGlobalGPUShapesVersion`, `tfx_GetGlobalColorRampsVersion`, `tfx_GetGPUGraphLookupsVersion`) each time it
  uploads. They change after a load, a free or a template customisation, then hold still. One that changes every frame
  means something is re-laying the libraries out each frame.
- Capture a frame, swap the order the libraries load in, and capture the same frame again with the same
  `tfx_SetStageSeed` and schedule. Swapping moves every library's place in the global buffers, so the two must match.
  Effects from one library drawing with another's shapes or colors means the backend is binding a library's own data
  (`tfx_GetParticlePropertiesBuffer`, `tfx_GetColorRampBitmap`, `tfx_GetLibraryGPUShapes`) instead of the global data.
- Free a library while another's effects keep playing, on a frame where the timer doesn't tick (pause first).
  Call `tfx_UpdateStage(stage, 0)` and upload the instances before drawing. The freed library's effects must be gone
  and the remaining effects must draw with their own library's data on that frame too. Drawing without the update shows
  them with the moved layout until the next real update.
- Then load it again while the others play and add its effects back. With every shape in one atlas texture (fixed
  bindings), it has to go back into its own layers, see "Freeing and reloading a library" in implementations/README.md.
  Its effects must draw with their own shapes and colors, and nothing else may flicker.
