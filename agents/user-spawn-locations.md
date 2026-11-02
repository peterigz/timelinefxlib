# User spawn locations: checks for agents

timelinefx.h is the reference (`tfx_AddSpawnLocation` and the functions after it). This file only holds the mistakes
that cost the most time.

**First:** read the log for lines starting `TimelineFX tip:`. If the host has no console, route them to its log with
`tfx_SetPrintCallback`.

## An effect that spawns nothing

An effect saved with "Spawn at User Locations" spawns nothing at all until it has a location, on the CPU and the GPU,
even though the stage counts it as running with live emitters. Every Run on GPU effect the editor saves is one of
these. Check before you rely on an effect's own position:

- `tfx_EffectTemplateSpawnsAtUserLocations(template)`, or `tfx_EffectSpawnsAtUserLocations(descriptor)` when walking a
  library.
- `tfx_PrintLibrary(library)` marks them `spawns at user spawn locations (tfx_AddSpawnLocation)`.
- Adding one to a stage prints a tip, and `tfx_DescribeStage` lists those with no location right now.

Then add the effect once and give it a location per source:

```c
tfxEffectID effect_id = tfx_AddEffectTemplateToStage(stage, effect_template);
//Each burst, ship or impact
tfxSpawnLocationID location = tfx_AddSpawnLocation(stage, effect_id, position, tfxSpawnLocationAdd_none);
```

Call these from the thread that calls `tfx_UpdateStage`. They're applied at the next update. On an effect saved without
"Spawn at User Locations", `tfx_AddSpawnLocation` does nothing and returns an id `tfx_SpawnLocationIsValid` rejects,
with a tip the first time.

## How long a location lives

- Single shot emitters fire once per location, so an effect made only of them plays once there. For something that
  follows a location for as long as it lives (a ship's engine glow or trail), use an effect whose emitters aren't single
  shot: `tfx_PrintLibrary` reads `each location plays until removed` for those.
- A finite effect's location removes itself once the effect has played out there. `tfx_GetEffectLifetime(template)` is
  how long that takes for one location, not for the whole effect. Pass `tfxSpawnLocationAdd_no_auto_remove` to keep it.
- The effect itself keeps running after its locations have gone, and `tfx_IsEffectExpired` stays false. For one-shot
  bursts, either:
  - add the effect once and give it a location per burst (best, one effect for all of them), or
  - add an effect per burst and pass `tfxSpawnLocationAdd_expire_effect`, which soft expires the effect once it has no
    locations left, so it ends with its last particle like any finite effect.
  A location kept with `tfxSpawnLocationAdd_no_auto_remove` never goes, so it never lets `_expire_effect` fire.
- `tfxSpawnLocationAdd_transient` removes the location after the next update, for one-off bursts where you don't keep
  the id.
- Ids are generational, so setters called with the id of a removed location do nothing.
- To work out how many locations play at once (for GPU ring sizes), divide `tfx_GetEffectLifetime` by how often you add
  one and add 1. See gpu-particles.md.

## Multipliers

`tfx_SetSpawnLocation*Multiplier` are applied as each particle spawns, so particles already alive keep their values.
The velocity multiplier scales base velocity, never the emitter's spawn impulse, so an impulse-driven burst barely
changes with it. That's true on the CPU and the GPU alike, and of `tfx_SetEffectVelocityMultiplier` too. Size and
life multipliers are a better way to make one burst look different from the next.
