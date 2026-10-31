# GPU particles: checks for agents

timelinefx.h is the reference (its "GPU particle simulation" notes and each function's comment). This file only holds
what to print, how to read it and the mistakes that cost the most time.

**First:** read the log for lines starting `TimelineFX tip:`. Each one prints once, names the effect or emitter, and
says what to change. If the host has no console, route them to its log with `tfx_SetPrintCallback`.

## Print the stage

After a second or two of updates, straight after a frame's `tfx_GetGPUParticleFrame`. The unfetched count is in ticks,
so printing between updating and fetching shows a false `Unfetched:` line whenever a frame runs more than 2:

```c
tfxU32 size = tfx_DescribeStage(stage, NULL, 0);
char *text = (char *)malloc(size);
tfx_DescribeStage(stage, text, size);
printf("%s", text);
free(text);
```

`tfx_GetStageDiagnostics` returns the same as a struct, and `tfx_GetGPUParticleStats` just the GPU numbers. To see
what a library contains (which emitters are GPU, ring sizes, which effects need spawn locations) call
`tfx_PrintLibrary(library)`.

## Reading it

| You see | It means | Do |
| --- | --- | --- |
| Rings allocated but `estimated_live_particles` 0, or a `No spawn locations right now` line | The effect isn't spawning. Usually it spawns at user spawn locations and has none | `tfx_AddSpawnLocation` (see user-spawn-locations.md) |
| `blocked_emitters` > 0, or a `Blocked:` line | The pool was too full to give the emitter a ring, so it spawns nothing | Size the pool from `tfx_SetTemplateGPUCapacityForLocations`, below |
| `dropped_spawns` > 0, or a `Dropping spawns:` line | The emitter's ring is too small for how many locations overlap | A larger location count in `tfx_SetTemplateGPUCapacityForLocations`, or `tfx_SetTemplateEmitterGPUCapacity` for one emitter |
| An `Unfetched:` line, or `unfetched_gpu_updates` climbing | The host updates the stage without calling `tfx_GetGPUParticleFrame` | Fetch and dispatch on every frame the stage updated |
| `slots_in_use` well above `estimated_live_particles` | Normal: a dead particle keeps its slot until every particle spawned in the same tick has died | Nothing. Compare `estimated_live_particles` with the CPU, never `slots_in_use` |
| A `Library ... tfxErrorCode_` line | The library loaded with problems, e.g. `file_version_out_of_date` | Re-save it in the editor |

## Sizing the pool

Every Run on GPU emitter instance gets a fixed ring, sized by default for one spawn location. For an effect that plays
at several locations at once:

```c
//How many locations play at once: one location's play time over how often you add one
tfxU32 overlap = (tfxU32)(tfx_GetEffectLifetime(effect_template) / milliseconds_between_locations) + 1;
stage_info.gpu_particle_capacity += tfx_SetTemplateGPUCapacityForLocations(effect_template, overlap);
//...the same for every template, then tfx_CreateStage(stage_info)
```

The return is already rounded the way the pool hands out rings, so the sum fits exactly. Don't add up
`tfx_GetEmitterGPUCapacity` yourself, which misses the rounding. Capacities assume about 60 updates a second. Fewer
updates means bigger ticks and drops.

## Template paths

Every `tfx_SetTemplate*` and `tfx_*TemplateEmitter` function takes `"Effect/Emitter"` or `tfx_GetDescriptorPath(emitter)`,
so walking `tfx_GetDescriptorChild` over a template's effect and passing each child's path works. A path that isn't in
the template prints a tip with a valid example, then asserts. `tfx_GetEmitterFromTemplate` checks a path without
asserting.

## Is it the effect or the integration?

`tfx_SetStageGPUSimulation(stage, false)`, then `tfx_ClearStage` and add the effects again (running effects keep the
simulation they started with). Still missing on the CPU: the effect isn't spawning, see above. Appears on the CPU: the
problem is in the host's upload, dispatches or instance ranges. Go through the GPU particles part of the porting
checklist in implementations/README.md.

## Run on GPU emitters that run on the CPU anyway

`tfx_PrintLibrary` prints `Run on GPU but runs on the CPU because ...` for these: other emitters spawn from them, or
their particles are relative and their effect doesn't spawn at user spawn locations. Every Run on GPU emitter also runs
on the CPU when the stage has no `gpu_particle_capacity` (there's a tip for that), while the stage records sprite data,
and while `tfx_SetStageGPUSimulation` has it off.
