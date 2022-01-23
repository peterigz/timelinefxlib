#include "timelinefx.h"
#include "Libraries/miniz.h"

namespace tfx {
	int FormatString(char* buf, size_t buf_size, const char* fmt, va_list args) {
		int w = vsnprintf(buf, buf_size, fmt, args);
		if (buf == NULL)
			return w;
		if (w == -1 || w >= (int)buf_size)
			w = (int)buf_size - 1;
		buf[w] = 0;
		return w;
	}

	void tfxText::Appendf(const char *format, ...) {
		va_list args;
		va_start(args, format);

		va_list args_copy;
		va_copy(args_copy, args);

		int len = FormatString(NULL, 0, format, args);         // FIXME-OPT: could do a first pass write attempt, likely successful on first pass.
		if (len <= 0)
		{
			va_end(args_copy);
			return;
		}

		const int write_off = (string.size() != 0) ? string.size() : 1;
		const int needed_sz = write_off + len;
		if (write_off + (unsigned int)len >= string.capacity)
		{
			int new_capacity = string.capacity * 2;
			string.reserve(needed_sz > new_capacity ? needed_sz : new_capacity);
		}

		string.resize(needed_sz);
		FormatString(&string[write_off - 1], (size_t)len + 1, format, args_copy);
		va_end(args_copy);

		va_end(args);
	}

	EffectEmitter::~EffectEmitter() {
		sub_effectors.free_all();
	}

	void EffectEmitter::SoftExpire() {
		flags |= tfxEmitterStateFlags_stop_spawning;
	}

	void EffectEmitter::Rotate(float r) {
		local.rotation += r;
	}
	void EffectEmitter::SetAngle(float a) {
		local.rotation = a;
	}
	void EffectEmitter::Scale(const tfxVec2& s) {
		local.scale = s;
	}
	void EffectEmitter::Scale(float x, float y) {
		local.scale.x = x;
		local.scale.y = y;
	}
	void EffectEmitter::Move(const tfxVec2& m) {
		local.position += m;
	}
	void EffectEmitter::Move(float x, float y) {
		local.position.x += x;
		local.position.y += y;
	}
	void EffectEmitter::Position(const tfxVec2& p) {
		local.position = p;
	}
	void EffectEmitter::Position(float x, float y) {
		local.position.x = x;
		local.position.y = y;
	}
	void EffectEmitter::UpdateMaxLife() {
		max_life = GetMaxLife(*this);
		GetGraphByType(tfxOvertime_red)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_green)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_blue)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_opacity)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_intensity)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_velocity)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_width)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_height)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_weight)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_spin)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_stretch)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_spin)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_frame_rate)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_motion_randomness)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_velocity_adjuster)->lookup.life = max_life;
		GetGraphByType(tfxOvertime_direction)->lookup.life = max_life;
	}

	EffectEmitter& EffectEmitter::AddEmitter(EffectEmitter &e) {
		e.type = EffectEmitterType::tfxEmitter;
		e.library = library;
		e.uid = ++library->uid;
		sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return sub_effectors.back();
	}

	EffectEmitter& EffectEmitter::AddEffect(EffectEmitter &e) {
		e.type = EffectEmitterType::tfxEffect;
		e.library = library;
		e.parent = this;
		e.uid = ++library->uid;
		sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return sub_effectors.back();
	}

	EffectEmitter& EffectEmitter::AddEffect() {
		EffectEmitter e;
		e.library = library;
		e.uid = ++library->uid;
		e.type = EffectEmitterType::tfxEffect;
		e.name = "New Effect";
		sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return sub_effectors.back();
	}

	EffectEmitter &EffectEmitter::AddEffector(EffectEmitterType type) {
		EffectEmitter e;
		//e.parent_effect = this;
		e.type = type;
		e.library = library;
		e.uid = ++library->uid;
		if(e.type == tfxEffect)
			e.name = "New Effect";
		else
			e.name = "New Emitter";
		sub_effectors.push_back(e);
		library->UpdateEffectPaths();
		ReIndex();
		return sub_effectors.back();
	}

	void EffectEmitter::Update() {
		captured = world;

		if (pm->lookup_mode == tfxPrecise) {
			current.frame = current.age;
		}
		else {
			current.frame = current.age / tfxLOOKUP_FREQUENCY;
		}

		if (type == EffectEmitterType::tfxEffect) {
			particle_count = 0;
			UpdateEffectState();
		}

		if (parent) {
			parent = parent->next_ptr;
			flags |= (parent->flags & tfxEmitterStateFlags_remove);
			UpdateEmitterState();
			TransformEffector(*parent, true, properties.flags & tfxEmitterPropertyFlags_relative_angle);

			if (flags & tfxEmitterStateFlags_no_tween_this_update) {
				captured = world;
			}

			if (!pm->disable_spawing && pm->FreeCapacity(properties.layer))
				SpawnParticles();
			parent->particle_count += particle_count;
		}
		else if (parent_particle) {
			flags |= parent_particle->flags & tfxParticleFlags_remove;
			if (parent_particle->next_ptr) {
				parent_particle = parent_particle->next_ptr;
				TransformEffector(*parent_particle, true, properties.flags & tfxEmitterPropertyFlags_relative_angle);
				world.position += properties.emitter_handle;

				if (flags & tfxEmitterStateFlags_no_tween_this_update) {
					captured = world;
				}

			}
			else {
				parent_particle = nullptr;
				flags |= tfxEmitterStateFlags_retain_matrix;
				local = world;
				flags |= tfxEmitterStateFlags_stop_spawning;
			}
		}
		else {
			if (!(flags & tfxEmitterStateFlags_retain_matrix)) {
				world = local;
				world.position += properties.emitter_handle;
				float s = sin(local.rotation);
				float c = cos(local.rotation);
				matrix.Set(c, s, -s, c);
			}

			if (flags & tfxEmitterStateFlags_no_tween_this_update) {
				captured = world;
			}
		}

		current.age += FRAME_LENGTH;

		if (properties.loop_length && current.age > properties.loop_length)
			current.age = 0;

		if (particle_count == 0) {
			timeout_counter++;
			if (parent && timeout_counter >= timeout)
				parent->active_children--;
		}
		else {
			timeout_counter = 0;
		}

		flags &= ~tfxEmitterStateFlags_no_tween_this_update;

	}

	void EffectEmitter::NoTweenNextUpdate() {
		tfxvec<EffectEmitter*> stack;
		stack.push_back(this);
		while (!stack.empty()) {
			auto &current = stack.pop_back();
			current->flags |= tfxEmitterStateFlags_no_tween_this_update;
			for (auto &sub : current->sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	void EffectEmitter::SpawnParticles() {
		if (current.single_shot_done || parent->flags & tfxEmitterStateFlags_stop_spawning)
			return;

		float qty;
		if (!(properties.flags & tfxEmitterPropertyFlags_single) && !(properties.flags & tfxEmitterPropertyFlags_one_shot)) {
			qty = current.amount;
			qty += random_generation.Range(current.amount_variation);
			qty *= lookup_callback(parent->library->global_graphs[parent->global].amount, current.frame);
			//qty *= parent->library->LookupPreciseNodeList(tfxGlobal_amount, parent->lookup_node_index, current.frame);
			//qty *= parent->library->LookupFastValueList(tfxGlobal_amount, parent->lookup_value_index, current.frame);
			qty *= UPDATE_TIME;
			qty += current.amount_remainder;
		}
		else {
			qty = (float)properties.spawn_amount;
		}

		float tween = 0.f;
		float interpolate = (float)(int)qty;
		float count = 0;

		while (qty >= 1.f) {
			if (!pm->FreeCapacity(properties.layer))
				break;
			tween = count / interpolate;
			count++;
			qty -= 1.f;

			Particle &p = pm->GrabParticle(properties.layer);
			p.flags = tfxParticleFlags_fresh;
			p.parent = this;
			p.next_ptr = &p;

			if (properties.flags & (tfxEmitterPropertyFlags_single | tfxEmitterPropertyFlags_one_shot))
				current.single_shot_done = true;

			//----Properties

			//Set base values-------------------------------

			//----Life
			p.max_age = current.life + random_generation.Range(current.life_variation);
			p.age = 0.f;

			//----Position
			if (properties.emission_type == EmissionType::tfxPoint) {
				if (properties.flags & tfxEmitterPropertyFlags_relative_position)
					p.local.position = -current.emitter_handle;
				else {
					if (properties.flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) {
						p.local.position = InterpolateVec2(tween, captured.position, world.position);
					}
					else {
						tfxVec2 rotvec = matrix.TransformVector(-current.emitter_handle);
						tfxVec2 spawn_position = InterpolateVec2(tween, captured.position, world.position) * world.scale;
						p.local.position = rotvec + spawn_position;
					}
				}
			}
			else if (properties.emission_type == EmissionType::tfxArea) {
				tfxVec2 position = tfxVec2(0.f, 0.f);

				if (properties.flags & tfxEmitterPropertyFlags_spawn_on_grid) {

					if (properties.flags & tfxEmitterPropertyFlags_fill_area) {
						if (!(properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
							current.grid_coords.x--;
							if (current.grid_coords.x < 0.f) {
								current.grid_coords.y--;
								current.grid_coords.x = properties.grid_points.x - 1;
								if (current.grid_coords.y < 0.f)
									current.grid_coords.y = properties.grid_points.y - 1;
							}
						}

						p.local.position = position + (current.grid_coords * current.grid_segment_size) + current.emitter_handle;

						if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
							current.grid_coords.x++;
							if (current.grid_coords.x == properties.grid_points.x) {
								current.grid_coords.y++;
								current.grid_coords.x = 0.f;
								if (current.grid_coords.y >= properties.grid_points.y)
									current.grid_coords.y = 0.f;
							}
						}
					}
					else {

						if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {

							current.grid_direction.x = 1;
							current.grid_direction.y = 0;
							if (current.grid_coords.x == properties.grid_points.x - 1 && current.grid_coords.y >= 0 && current.grid_coords.y < properties.grid_points.y - 1) {
								current.grid_direction.x = 0;
								current.grid_direction.y = 1;
							}
							else if (current.grid_coords.x > 0 && current.grid_coords.x < properties.grid_points.x && current.grid_coords.y == properties.grid_points.y - 1) {
								current.grid_direction.x = -1;
								current.grid_direction.y = 0;
							}
							else if (current.grid_coords.x == 0 && current.grid_coords.y > 0 && current.grid_coords.y < properties.grid_points.y) {
								current.grid_direction.x = 0;
								current.grid_direction.y = -1;
							}

						}
						else {

							current.grid_direction.x = -1;
							current.grid_direction.y = 0;
							if (current.grid_coords.x == properties.grid_points.x - 1 && current.grid_coords.y > 0 && current.grid_coords.y < properties.grid_points.y) {
								current.grid_direction.x = 0;
								current.grid_direction.y = -1;
							}
							else if (current.grid_coords.x >= 0 && current.grid_coords.x < properties.grid_points.x - 1 && current.grid_coords.y == properties.grid_points.y - 1) {
								current.grid_direction.x = 1;
								current.grid_direction.y = 0;
							}
							else if (current.grid_coords.x == 0 && current.grid_coords.y >= 0 && current.grid_coords.y < properties.grid_points.y - 1) {
								current.grid_direction.x = 0;
								current.grid_direction.y = 1;
							}

						}

						current.grid_coords += current.grid_direction; 
						tfxBound(current.grid_coords, properties.grid_points);
						p.local.position = position + (current.grid_coords * current.grid_segment_size) + current.emitter_handle;
					}
				}
				else {
					if (properties.flags & tfxEmitterPropertyFlags_fill_area) {
						position.x = random_generation.Range(current.emitter_size.x);
						position.y = random_generation.Range(current.emitter_size.y);
					}
					else {
						//Spawn on one of 4 edges of the area
						unsigned int side = random_generation.RangeUInt(4);
						if (side == 0) {
							//left side
							position.x = 0.f;
							position.y = random_generation.Range(current.emitter_size.y);
						}
						else if (side == 1) {
							//right side
							position.x = current.emitter_size.x;
							position.y = random_generation.Range(current.emitter_size.y);
						}
						else if (side == 2) {
							//top side
							position.x = random_generation.Range(current.emitter_size.x);
							position.y = 0.f;
						}
						else if (side == 3) {
							//bottom side
							position.x = random_generation.Range(current.emitter_size.x);
							position.y = current.emitter_size.y;
						}
					}

					p.local.position = position + current.emitter_handle;
				}

				//----TForm and Emission
				if (!(properties.flags & tfxEmitterPropertyFlags_relative_position)) {
					p.local.position = matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));
					p.local.position = world.position + p.local.position * world.scale;
				}

			}
			else if (properties.emission_type == EmissionType::tfxEllipse) {
				tfxVec2 emitter_size = (current.emitter_size * .5f);
				tfxVec2 position = tfxVec2(0.f, 0.f);

				if (properties.flags & tfxEmitterPropertyFlags_spawn_on_grid && !(properties.flags & tfxEmitterPropertyFlags_fill_area)) {

					current.grid_coords.y = 0.f;

					if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						current.grid_coords.x--;
						if (current.grid_coords.x < 0.f) {
							current.grid_coords.x = properties.grid_points.x - 1;
						}
					}

					float th = current.grid_coords.x * current.grid_segment_size.x + current.arc_offset;
					p.local.position = tfxVec2(std::cosf(th) * emitter_size.x + current.emitter_handle.x + emitter_size.x,
						-std::sinf(th) * emitter_size.y + current.emitter_handle.y + emitter_size.y);

					if (!(properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						current.grid_coords.x++;
						if (current.grid_coords.x >= properties.grid_points.x) {
							current.grid_coords.x = 0.f;
						}
					}

				}
				else if(!(properties.flags & tfxEmitterPropertyFlags_fill_area)) {
					float th = random_generation.Range(current.arc_size) + current.arc_offset;

					p.local.position = tfxVec2(std::cosf(th) * emitter_size.x + current.emitter_handle.x + emitter_size.x,
						-std::sinf(th) * emitter_size.y + current.emitter_handle.y + emitter_size.y);

				}
				else {
					p.local.position.x = random_generation.Range(-emitter_size.x, emitter_size.x);
					p.local.position.y = random_generation.Range(-emitter_size.y, emitter_size.y);

					while ((std::pow(p.local.position.x, 2) / std::pow(emitter_size.x, 2)) + (std::pow(p.local.position.y, 2) / std::pow(emitter_size.y, 2)) > 1) {
						p.local.position.x = random_generation.Range(-emitter_size.x, emitter_size.x);
						p.local.position.y = random_generation.Range(-emitter_size.y, emitter_size.y);
					}
				}

				//----TForm and Emission
				if (!(properties.flags & tfxEmitterPropertyFlags_relative_position)) {
					p.local.position = matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));
					p.local.position = world.position + p.local.position * world.scale;
				}

			}
			else if (properties.emission_type == EmissionType::tfxLine) {
				if (properties.flags & tfxEmitterPropertyFlags_spawn_on_grid) {

					current.grid_coords.x = 0.f;

					if (!(properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise)) {
						current.grid_coords.y--;
						if (current.grid_coords.y < 0.f) {
							current.grid_coords.y = properties.grid_points.x - 1;
						}
					}

					p.local.position = tfxVec2(current.grid_coords * -current.grid_segment_size);
					p.distance_travelled = std::fabsf(p.local.position.y);
					p.local.position += current.emitter_handle;

					if (properties.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) {
						current.grid_coords.y++;
						if (current.grid_coords.y >= properties.grid_points.x) {
							current.grid_coords.y = 0.f;
						}
					}

				}
				else {
					p.local.position.x = 0.f;
					p.local.position.y = random_generation.Range(-current.emitter_size.y, 0.f);
					p.distance_travelled = std::fabsf(p.local.position.y);

					p.local.position += current.emitter_handle;

				}

				//----TForm and Emission
				if (!(properties.flags & tfxEmitterPropertyFlags_relative_position) && !(properties.flags & tfxEmitterPropertyFlags_edge_traversal)) {
					p.local.position = matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));
					p.local.position = world.position + p.local.position * world.scale;
				}

			}

			//----Weight
			if (current.weight) {
				p.base.weight = current.weight * library->overtime_graphs[overtime].weight.GetFirstValue();
				if (current.weight_variation > 0) {
					p.base.weight += random_generation.Range(-current.weight_variation, current.weight_variation) * library->overtime_graphs[overtime].weight.GetFirstValue();
				}
			}
			else {
				p.base.weight = 0;
			}
			p.weight_acceleration = p.base.weight * library->overtime_graphs[overtime].weight.GetFirstValue() * UPDATE_TIME;

			//----Velocity
			p.velocity_scale = library->overtime_graphs[overtime].velocity.GetFirstValue() * current.velocity_adjuster;
			p.base.velocity = current.velocity + random_generation.Range(-current.velocity_variation, current.velocity_variation);

			//----Size
			if (!(properties.flags & tfxEmitterPropertyFlags_base_uniform_size)) {
				p.base.random_size.x = random_generation.Range(current.size_variation.x);
				p.base.random_size.y = random_generation.Range(current.size_variation.y);
				p.base.height = p.base.random_size.y + current.size.y;
				p.base.size.x = (p.base.random_size.x + current.size.x) / properties.image->image_size.x;
				p.base.size.y = p.base.height / properties.image->image_size.y;

				p.local.scale.x = p.base.size.x * library->overtime_graphs[overtime].width.GetFirstValue();

				if (library->overtime_graphs[overtime].stretch.GetFirstValue()) {
					float velocity = std::fabsf(p.velocity_scale * p.base.velocity) * UPDATE_TIME;
					velocity += p.weight_acceleration * UPDATE_TIME;
					p.local.scale.y = (library->overtime_graphs[overtime].height.GetFirstValue() * parent->current.size.y * (p.base.height + (velocity * library->overtime_graphs[overtime].stretch.GetFirstValue() * parent->current.stretch))) / properties.image->image_size.y;
				}
				else {
					if (properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
						p.local.scale.y = p.base.size.y * library->overtime_graphs[overtime].width.GetFirstValue();
					}
					else {
						p.local.scale.y = p.base.size.y * library->overtime_graphs[overtime].height.GetFirstValue();
					}
				}
			}
			else {
				p.base.random_size.x = random_generation.Range(current.size_variation.x);
				p.base.random_size.y = p.base.random_size.x;
				p.base.height = p.base.random_size.y + current.size.y;
				p.base.size.x = (p.base.random_size.x + current.size.x) / properties.image->image_size.x;
				p.base.size.y = p.base.height / properties.image->image_size.y;

				p.local.scale.x = p.base.size.x * library->overtime_graphs[overtime].width.GetFirstValue();

				if (library->overtime_graphs[overtime].stretch.GetFirstValue()) {
					float velocity = std::fabsf(p.velocity_scale * p.base.velocity) * UPDATE_TIME;
					velocity += p.weight_acceleration * UPDATE_TIME;
					p.local.scale.y = (library->overtime_graphs[overtime].width.GetFirstValue() * parent->current.size.y * (p.base.height + (velocity * library->overtime_graphs[overtime].stretch.GetFirstValue() * parent->current.stretch))) / properties.image->image_size.y;
				}
				else {
					p.local.scale.y = p.local.scale.x;
				}
			}

			//----Spin
			p.base.spin = random_generation.Range(-current.spin_variation, std::abs(current.spin_variation)) + current.spin;

			switch (properties.angle_setting) {
			case AngleSetting::tfxRandom:
				p.local.rotation = random_generation.Range(properties.angle_offset);
				break;
			case AngleSetting::tfxSpecify:
				p.local.rotation = properties.angle_offset;
				break;
			default:
				break;
			}

			//----Splatter
			if (current.splatter) {
				float splattertemp = current.splatter;
				float splatx = random_generation.Range(-current.splatter, current.splatter);
				float splaty = random_generation.Range(-current.splatter, current.splatter);

				while (GetDistance(0, 0, splatx, splaty) >= splattertemp && splattertemp > 0) {
					splatx = random_generation.Range(-current.splatter, current.splatter);
					splaty = random_generation.Range(-current.splatter, current.splatter);
				}

				if (!(properties.flags & tfxEmitterPropertyFlags_relative_position)) {
					p.local.position.x += splatx * world.scale.x;
					p.local.position.y += splaty * world.scale.y;
				}
				else {
					p.local.position.x += splatx;
					p.local.position.y += splaty;
				}
			}

			float direction = 0;

			if (properties.angle_setting == AngleSetting::tfxAlign && properties.flags & tfxEmitterPropertyFlags_edge_traversal)
				p.world.rotation = p.local.rotation = direction + properties.angle_offset;

			TransformParticlePrevious(p, *this);
			p.captured = p.world;
			TransformParticle(p, *this);

			//----Motion randomness
			p.motion_randomness = current.motion_randomness;
			float mr = tfxRadians(p.motion_randomness * library->overtime_graphs[overtime].motion_randomness.GetFirstValue());
			float motion_randomness_direction = tfxRadians(22.5f * random_generation.Range(-mr, mr));
			p.motion_randomness_speed = 30.f * random_generation.Range(-mr, mr);

			if (!(properties.flags & tfxEmitterPropertyFlags_edge_traversal) || properties.emission_type != EmissionType::tfxLine) {
				direction = p.emission_angle = GetEmissionDirection(p) + motion_randomness_direction + library->overtime_graphs[overtime].direction.GetFirstValue();
			}

			//----Normalize Velocity to direction
			tfxVec2 velocity_normal;
			velocity_normal.x = std::sinf(direction);
			velocity_normal.y = -std::cosf(direction);

			//p.velocity = p.velocity_normal * p.base.velocity * p.velocity_scale * UPDATE_TIME;
			bool line = properties.flags & tfxEmitterPropertyFlags_edge_traversal && properties.emission_type == EmissionType::tfxLine;

			if (properties.angle_setting == AngleSetting::tfxAlign && !line) {
				p.world.rotation = p.local.rotation = GetVectorAngle(velocity_normal.x, velocity_normal.y) + properties.angle_offset;
				if (properties.flags & tfxEmitterPropertyFlags_relative_angle)
					p.world.rotation += world.rotation;
				p.captured.rotation = p.world.rotation;
			}

			//----Handle
			/*if (properties.flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
				p.handle = tfxVec2(0.5f, 0.5f);
			}
			else {
				p.handle = properties.image_handle;
			}*/

			//----Image
			//p.image = properties.image;
			if (properties.flags & tfxEmitterPropertyFlags_random_start_frame && properties.image->animation_frames > 1) {
				p.image_frame = random_generation.Range(properties.image->animation_frames);
			}
			else {
				p.image_frame = properties.start_frame;
			}
			p.image_frame_rate = library->overtime_graphs[overtime].frame_rate.GetFirstValue();

			//----Color
			p.color.a = unsigned char(255.f * library->overtime_graphs[overtime].opacity.GetFirstValue() * parent->current.color.a);
			p.intensity = library->overtime_graphs[overtime].opacity.GetFirstValue();
			if (properties.flags & tfxEmitterPropertyFlags_random_color) {
				float age = random_generation.Range(p.max_age);
				p.color.r = unsigned char(255.f * lookup_overtime_callback(library->overtime_graphs[overtime].red, age, p.max_age));
				p.color.g = unsigned char(255.f * lookup_overtime_callback(library->overtime_graphs[overtime].green, age, p.max_age));
				p.color.b = unsigned char(255.f * lookup_overtime_callback(library->overtime_graphs[overtime].blue, age, p.max_age));
			}
			else {
				p.color.r = unsigned char(255.f * library->overtime_graphs[overtime].red.GetFirstValue());
				p.color.g = unsigned char(255.f * library->overtime_graphs[overtime].green.GetFirstValue());
				p.color.b = unsigned char(255.f * library->overtime_graphs[overtime].blue.GetFirstValue());
			}

			if (sub_effectors.size()) {
				for (auto &e : sub_effectors) {
					if (!pm->FreeEffectCapacity())
						break;
					e.parent = nullptr;
					e.parent_particle = &p;
					e.particle_count = 0;
					pm->AddEffect(e, !pm->current_ebuff);
				}
			}

			if (particle_onspawn_callback)
				particle_onspawn_callback(p);

			particle_count++;
		}

		current.amount_remainder = qty;
	}

	float EffectEmitter::GetEmissionDirection(Particle& p) {

		//----Emission
		float range = current.emission_angle_variation *.5f;
		float direction = 0;

		if (properties.emission_type == EmissionType::tfxPoint)
			return direction + current.emission_angle + random_generation.Range(-range, range);

		tfxVec2 tmp_position;
		if (p.local.position.x == 0 && p.local.position.y == 0)
			tmp_position = current.emitter_size;
		else
			tmp_position = p.local.position;

		if (properties.emission_direction == EmissionDirection::tfxOutwards) {

			tfxVec2 to_handle;

			if (properties.flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = (tmp_position);
			else
				to_handle = (p.world.position - world.position);

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (properties.emission_direction == EmissionDirection::tfxInwards) {

			tfxVec2 to_handle;

			if (properties.flags & tfxEmitterPropertyFlags_relative_position)
				to_handle = (-tmp_position);
			else
				to_handle = (world.position - p.world.position);

			direction = GetVectorAngle(to_handle.x, to_handle.y);

		}
		else if (properties.emission_direction == EmissionDirection::tfxBothways) {

			if (current.emission_alternator) {

				tfxVec2 to_handle;

				if (properties.flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = (tmp_position);
				else
					to_handle = (p.world.position - world.position);

				direction = GetVectorAngle(to_handle.x, to_handle.y);

			}
			else {

				tfxVec2 to_handle;

				if (properties.flags & tfxEmitterPropertyFlags_relative_position)
					to_handle = (-tmp_position);
				else
					to_handle = (world.position - p.world.position);

				direction = GetVectorAngle(to_handle.x, to_handle.y);

			}

			current.emission_alternator = !current.emission_alternator;
		}

		if (std::isnan(direction))
			direction = 0.f;
		return direction + current.emission_angle + random_generation.Range(-range, range);
	}

	void EffectEmitter::UpdateEmitterState() {
		//todo: tidy this up, and not all of this needs to be called every frame
		if (parent == nullptr) {
			//is this even valid now? Every emitter should have an effect as a parent
			current.life = lookup_callback(library->base_graphs[base].life, current.frame);
			current.life_variation = lookup_callback(library->variation_graphs[variation].life, current.frame);
			if (!(properties.flags & tfxEmitterPropertyFlags_base_uniform_size)) {
				current.size.x = lookup_callback(library->base_graphs[base].width, current.frame);
				current.size.y = lookup_callback(library->base_graphs[base].height, current.frame);
			}
			else {
				current.size.x = lookup_callback(library->base_graphs[base].width, current.frame);
				current.size.y = current.size.x;
			}
			current.size_variation.x = lookup_callback(library->variation_graphs[variation].width, current.frame);
			current.size_variation.y = lookup_callback(library->variation_graphs[variation].height, current.frame);
			current.velocity = lookup_callback(library->base_graphs[base].velocity, current.frame);
			current.velocity_variation = lookup_callback(library->variation_graphs[variation].velocity, current.frame);
			current.velocity_adjuster = lookup_callback(library->overtime_graphs[overtime].velocity_adjuster, current.frame);
			current.spin = lookup_callback(library->base_graphs[base].spin, current.frame);
			current.spin_variation = lookup_callback(library->variation_graphs[variation].spin, current.frame);
			local.rotation = lookup_callback(library->property_graphs[property].emitter_angle, current.frame);
			current.emission_angle = lookup_callback(library->property_graphs[property].emission_angle, current.frame);
			current.emission_angle_variation = lookup_callback(library->property_graphs[property].emission_range, current.frame);
			current.color.r = library->overtime_graphs[overtime].red.GetFirstValue();
			current.color.g = library->overtime_graphs[overtime].green.GetFirstValue();
			current.color.b = library->overtime_graphs[overtime].blue.GetFirstValue();
			current.color.a = library->overtime_graphs[overtime].opacity.GetFirstValue();
			current.splatter = lookup_callback(library->property_graphs[property].splatter, current.frame);
			current.emitter_size.x = lookup_callback(library->property_graphs[property].emitter_width, current.frame);
			current.weight = lookup_callback(library->base_graphs[base].weight, current.frame);
			current.weight_variation = lookup_callback(library->variation_graphs[variation].weight, current.frame);
			current.motion_randomness = lookup_callback(library->variation_graphs[variation].motion_randomness, current.frame);

			if (properties.emission_type == EmissionType::tfxArea || properties.emission_type == EmissionType::tfxEllipse)
				current.emitter_size.y = lookup_callback(library->property_graphs[property].emitter_width, current.frame);
			else
				current.emitter_size.y = 0.f;

			if (properties.emission_type == EmissionType::tfxEllipse) {
				current.arc_size = lookup_callback(library->property_graphs[property].arc_size, current.frame);
				current.arc_offset = lookup_callback(library->property_graphs[property].arc_offset, current.frame);
			}

		}
		else {
			EffectEmitter &e = *parent;
			current.amount = lookup_callback(library->base_graphs[base].amount, current.frame);
			current.amount_variation = lookup_callback(library->variation_graphs[base].amount, current.frame);
			current.life = lookup_callback(library->base_graphs[base].life, current.frame) * e.current.life;
			current.life_variation = lookup_callback(library->variation_graphs[variation].life, current.frame) * e.current.life;
			if (!(properties.flags & tfxEmitterPropertyFlags_base_uniform_size)) {
				current.size.x = lookup_callback(library->base_graphs[base].width, current.frame) * e.current.size.x;
				current.size.y = lookup_callback(library->base_graphs[base].height, current.frame) * e.current.size.y;
			}
			else {
				current.size.x = lookup_callback(library->base_graphs[base].width, current.frame);
				if (e.properties.flags & tfxEmitterPropertyFlags_global_uniform_size)
					current.size.y = current.size.x * e.current.size.x;
				else
					current.size.y = current.size.x * e.current.size.y;
				current.size.x *= e.current.size.x;
			}
			current.size_variation.x = lookup_callback(library->variation_graphs[variation].width, current.frame) * e.current.size.x;
			current.size_variation.y = lookup_callback(library->variation_graphs[variation].height, current.frame) * e.current.size.y;
			current.velocity = lookup_callback(library->base_graphs[base].velocity, current.frame) * e.current.velocity;
			current.velocity_variation = lookup_callback(library->variation_graphs[variation].velocity, current.frame) * e.current.velocity;
			current.velocity_adjuster = lookup_callback(library->overtime_graphs[overtime].velocity_adjuster, current.frame);
			current.spin = lookup_callback(library->base_graphs[base].spin, current.frame) * e.current.spin;
			current.spin_variation = lookup_callback(library->variation_graphs[variation].spin, current.frame) * e.current.spin;
			local.rotation = lookup_callback(library->property_graphs[property].emitter_angle, current.frame);
			current.emission_angle = lookup_callback(library->property_graphs[property].emission_angle, current.frame);
			current.emission_angle_variation = lookup_callback(library->property_graphs[property].emission_range, current.frame);
			current.color.r = library->overtime_graphs[overtime].red.GetFirstValue();
			current.color.g = library->overtime_graphs[overtime].green.GetFirstValue();
			current.color.b = library->overtime_graphs[overtime].blue.GetFirstValue();
			current.color.a = e.current.color.a;
			current.splatter = lookup_callback(library->property_graphs[property].splatter, current.frame) * e.current.splatter;
			current.emitter_size.y = lookup_callback(library->property_graphs[property].emitter_height, current.frame);
			current.weight = lookup_callback(library->base_graphs[base].weight, current.frame) * e.current.weight;
			current.weight_variation = lookup_callback(library->variation_graphs[variation].weight, current.frame) * e.current.weight;
			current.motion_randomness = lookup_callback(library->variation_graphs[variation].motion_randomness, current.frame);
			current.stretch = e.current.stretch;
			local.scale = e.local.scale;

			//----Handle
			if (properties.flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
				current.image_handle = tfxVec2(0.5f, 0.5f);
			}
			else {
				current.image_handle = properties.image_handle;
			}

			if (properties.emission_type == EmissionType::tfxArea || properties.emission_type == EmissionType::tfxEllipse) {
				current.emitter_size.x = lookup_callback(library->property_graphs[property].emitter_width, current.frame);
			}
			else
				current.emitter_size.x = 0.f;

			if (properties.emission_type == EmissionType::tfxEllipse) {
				current.arc_size = lookup_callback(library->property_graphs[property].arc_size, current.frame);
				current.arc_offset = lookup_callback(library->property_graphs[property].arc_offset, current.frame);
			}

		}

		if (properties.flags & tfxEmitterPropertyFlags_emitter_handle_auto_center && properties.emission_type != EmissionType::tfxPoint) {
			current.emitter_handle = current.emitter_size * -0.5f;
		}
		else if (!(properties.flags & tfxEmitterPropertyFlags_emitter_handle_auto_center)) {
			current.emitter_handle = properties.emitter_handle;
		}
		else {
			current.emitter_handle.x = current.emitter_handle.y = 0;
		}

		if (properties.flags & tfxEmitterPropertyFlags_emitter_handle_auto_center && properties.emission_type == EmissionType::tfxLine) {
			current.emitter_handle = current.emitter_size * 0.5f;
		}

		if (properties.flags & tfxEmitterPropertyFlags_spawn_on_grid) {
			if (properties.emission_type == EmissionType::tfxArea) {
				if (properties.grid_points.x > 1)
					current.grid_segment_size.x = current.emitter_size.x / (properties.grid_points.x - 1);
				if (properties.grid_points.y > 1)
					current.grid_segment_size.y = current.emitter_size.y / (properties.grid_points.y - 1);
			}
			else if (properties.emission_type == EmissionType::tfxEllipse) {
				if (properties.grid_points.x > 0)
					current.grid_segment_size.x = current.arc_size / (properties.grid_points.x);
			}
			else if (properties.emission_type == EmissionType::tfxLine) {
				if (properties.grid_points.x > 1)
					current.grid_segment_size.y = current.emitter_size.y / (properties.grid_points.x - 1);
			}
		}

		if (update_callback)
			update_callback(*this);
	}

	void EffectEmitter::UpdateEffectState() {
		//If this effect is a sub effect then the graph index will reference the global graphs for the root parent effect
		current.life = lookup_callback(library->global_graphs[global].life, current.frame);
		current.amount = lookup_callback(library->global_graphs[global].amount, current.frame);
		if (!(properties.flags & tfxEmitterPropertyFlags_global_uniform_size)) {
			current.size.x = lookup_callback(library->global_graphs[global].width, current.frame);
			current.size.y = lookup_callback(library->global_graphs[global].height, current.frame);
		}
		else {
			current.size.x = lookup_callback(library->global_graphs[global].width, current.frame);
			current.size.y = current.size.x;
		}
		current.velocity = lookup_callback(library->global_graphs[global].velocity, current.frame);
		current.spin = lookup_callback(library->global_graphs[global].spin, current.frame);
		current.color.a = lookup_callback(library->global_graphs[global].opacity, current.frame);
		current.splatter = lookup_callback(library->global_graphs[global].splatter, current.frame);
		//We don't want to scale twice when the sub effect is transformed, so the values here are set to 1. That means that the root effect will only control the global scale.
		if (!parent_particle) {
			local.scale.x = lookup_callback(library->global_graphs[global].overal_scale, current.frame);
			local.scale.y = local.scale.x;
			local.rotation = lookup_callback(library->global_graphs[global].effect_angle, current.frame);
		}
		else {
			local.scale.x = 1.f;
			local.scale.y = 1.f;
			local.rotation = 0.f;
		}
		current.stretch = lookup_callback(library->global_graphs[global].stretch, current.frame);
		current.weight = lookup_callback(library->global_graphs[global].weight, current.frame);

		if (update_callback)
			update_callback(*this);
	}

	void ReloadBaseValues(Particle &p, EffectEmitter &e) {
		//----Life
		//std::uniform_real_distribution<float> random_life(0, e.current.life_variation);
		//p.max_age = e.current.life + random_life(random_generation.engine);

		//----Velocity
		//std::uniform_real_distribution<float> random_velocity;
		//if (e.current.velocity_variation > 0)
			//random_velocity = std::uniform_real_distribution<float>(0, e.current.velocity_variation);
		//else
			//random_velocity = std::uniform_real_distribution<float>(e.current.velocity_variation, 0);
		//p.velocity_scale = e.library->overtime_graphs[e.overtime].velocity.GetFirstValue() * e.current.velocity_adjuster;
		//p.base.velocity = e.current.velocity + random_velocity(random_generation.engine);

		//----Size
		if (!(e.properties.flags & tfxEmitterPropertyFlags_base_uniform_size)) {
			//std::uniform_real_distribution<float> random_width(0, e.current.size_variation.x);
			//std::uniform_real_distribution<float> random_height(0, e.current.size_variation.y);

			p.base.height = p.base.random_size.y + e.current.size.y;
			p.base.size.x = (p.base.random_size.x + e.current.size.x) / e.properties.image->image_size.x;
			p.base.size.y = p.base.height / e.properties.image->image_size.y;

			//p.local.scale.x = p.base.size.x * e.library->overtime_graphs[e.overtime].width.GetFirstValue();

			if (e.library->overtime_graphs[e.overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(p.velocity_scale * p.base.velocity);
				//p.local.scale.y = (e.library->overtime_graphs[e.overtime].height.GetFirstValue() * e.parent->current.size.y * (p.base.height + (velocity * e.library->overtime_graphs[e.overtime].stretch.GetFirstValue() * e.parent->current.stretch))) / e.properties.image.image_size.y;
			}
			else {
				if (e.properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
					//p.local.scale.y = p.base.size.y * e.library->overtime_graphs[e.overtime].width.GetFirstValue();
				}
				else {
					//p.local.scale.y = p.base.size.y * e.library->overtime_graphs[e.overtime].height.GetFirstValue();
				}
			}
		}
		else {
			//std::uniform_real_distribution<float> random_width(0, e.current.size_variation.x);

			p.base.height = p.base.random_size.y + e.current.size.x;
			p.base.size.x = (p.base.random_size.x + e.current.size.x) / e.properties.image->image_size.x;
			p.base.size.y = p.base.height / e.properties.image->image_size.y;

			//p.local.scale.x = p.base.size.x * e.library->overtime_graphs[e.overtime].width.GetFirstValue();

			if (e.library->overtime_graphs[e.overtime].stretch.GetFirstValue()) {
				float velocity = std::fabsf(p.velocity_scale * p.base.velocity);
				//p.local.scale.y = (e.library->overtime_graphs[e.overtime].height.GetFirstValue() * e.parent->current.size.y * (p.base.height + (velocity * e.library->overtime_graphs[e.overtime].stretch.GetFirstValue() * e.parent->current.stretch))) / e.properties.image.image_size.y;
			}
			else {
				//p.local.scale.y = p.local.scale.x;
			}
		}

		//----Spin
		p.base.spin = random_generation.Range(-e.current.spin_variation, std::abs(e.current.spin_variation)) + e.current.spin;

		//std::uniform_real_distribution<float> random_angle(0, e.properties.angle_offset);
		//switch (e.properties.angle_setting) {
		//case AngleSetting::kRandom:
			//p.local.rotation = random_angle(random_generation.engine);
			//break;
		//case AngleSetting::kSpecify:
			//p.local.rotation = e.properties.angle_offset;
			//break;
		//default:
			//break;
		//}

		//----Weight
		if (e.current.weight) {
			if (e.current.weight_variation > 0) {
				p.base.weight = (random_generation.Range(e.current.weight_variation) + e.current.weight) * e.library->overtime_graphs[e.overtime].weight.GetFirstValue();
			}
			else {
				p.base.weight = (random_generation.Range(e.current.weight_variation, 0) + e.current.weight) * e.library->overtime_graphs[e.overtime].weight.GetFirstValue();
			}
			//p.weight_acceleration = 0;
		}
		else {
			//p.weight_acceleration = 0;
			p.base.weight = 0;
		}

		//TransformParticlePrevious(p, e);
		//p.captured = p.world;
		//TransformParticle(p, e);

		//----Motion randomness
		//p.motion_randomness = e.current.motion_randomness;
		//float mr = tfxRadians(p.motion_randomness * e.library->overtime_graphs[e.overtime].motion_randomness.GetFirstValue());
		//std::uniform_real_distribution<float> random_motion(-mr, mr);
		//p.motion_randomness_direction = tfxRadians(22.5f * random_motion(random_generation.engine));
		//p.motion_randomness_speed = 30.f * random_motion(random_generation.engine);
		//p.motion_tracker = 0;

		//if (!e.properties.edge_traversal || e.properties.emission_type != EmissionType::kLine) {
			//p.direction = p.emission_angle = e.GetEmissionDirection(p) + p.motion_randomness_direction;
		//}

		//----Normalize Velocity to direction
		//p.velocity_normal.x = std::sinf(p.direction);
		//p.velocity_normal.y = -std::cosf(p.direction);

		//p.velocity = p.velocity_normal * p.base.velocity * p.velocity_scale * e.timer->UpdateTime();
		//bool line = e.properties.edge_traversal && e.properties.emission_type == EmissionType::kLine;

		//if (e.properties.angle_setting == AngleSetting::kAlign && !line) {
			//p.world.rotation = p.local.rotation = GetVectorAngle(p.velocity_normal.x, p.velocity_normal.y) + e.properties.angle_offset + e.world.rotation;
			//p.captured.rotation = p.world.rotation;
		//}

		//----Handle
		//if (e.properties.flags & tfxEmitterPropertyFlags_image_handle_auto_center) {
			//p.handle = tfxVec2(0.5f, 0.5f);
		//}
		//else {
			//p.handle = e.properties.image_handle;
		//}

		//----Image
		//p.image = properties.image;
		//if (e.properties.image.random_start_frame && e.properties.image.animation_frames > 1) {
			//std::uniform_real_distribution<float> random_start_frame(0.f, e.properties.image.animation_frames - 1);
			//p.image_frame = random_start_frame(random_generation.engine);
		//}
		//else {
			//p.image_frame = e.properties.image.start_frame;
		//}
		//p.image_frame_rate = e.library->overtime_graphs[e.overtime].frame_rate.GetFirstValue();

		//----Color
		//p.color.a = unsigned char(255.f * e.library->overtime_graphs[e.overtime].opacity.GetFirstValue());
		//p.intensity = e.library->overtime_graphs[e.overtime].opacity.GetFirstValue();
		//if (e.properties.random_color) {
			//std::uniform_real_distribution<float> random_color(0.f, p.max_age);
			//float age = random_color(random_generation.engine);
			//p.color.r = unsigned char(255.f * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].red, age, p.max_age));
			//p.color.g = unsigned char(255.f * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].green, age, p.max_age));
			//p.color.b = unsigned char(255.f * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].blue, age, p.max_age));
		//}
		//else {
			//p.color.r = unsigned char(255.f * e.library->overtime_graphs[e.overtime].red.GetFirstValue());
			//p.color.g = unsigned char(255.f * e.library->overtime_graphs[e.overtime].green.GetFirstValue());
			//p.color.b = unsigned char(255.f * e.library->overtime_graphs[e.overtime].blue.GetFirstValue());
		//}

	}

	bool ControlParticle(Particle &p, EffectEmitter &e) {
		if (e.pm->update_base_values)
			ReloadBaseValues(p, e);

		p.age += FRAME_LENGTH;

		//-------------------------------------------------------
		//Controll what the particle does over the course of
		//it's lifetime
		//-------------------------------------------------------

		//Before we do anything, see if the particle should be remove/end of life
		p.flags |= e.flags & tfxParticleFlags_remove;
		if (p.flags & tfxParticleFlags_remove)
			return false;

		if (p.age >= p.max_age) {
			if (e.properties.flags & tfxEmitterPropertyFlags_single && !(e.properties.flags & tfxEmitterPropertyFlags_one_shot) && !e.pm->disable_spawing)
				p.age = 0.f;
			else {
				return false;
			}
		}
		float direction = 0;
		if (e.properties.emission_type != tfxLine && !(e.properties.flags & tfxEmitterPropertyFlags_edge_traversal)) {
			//----Motion randomness
			float mr = p.motion_randomness * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].motion_randomness, p.age, p.max_age);
			if (!(unsigned int)std::fmodf(p.age, FRAME_LENGTH * 2.f)) {
				p.motion_randomness_speed += 30.f * random_generation.Range(-mr, mr);
				p.emission_angle += tfxRadians(22.5f * random_generation.Range(-mr, mr));
			}
			direction = p.emission_angle + lookup_overtime_callback(e.library->overtime_graphs[e.overtime].direction, p.age, p.max_age);
		}

		//----Weight Changes
		p.weight_acceleration += p.base.weight * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].weight, p.age, p.max_age) * UPDATE_TIME;

		//----Velocity Changes
		tfxVec2 velocity_normal;
		velocity_normal.x = std::sinf(direction);
		velocity_normal.y = -std::cosf(direction);
		p.velocity_scale = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].velocity, p.age, p.max_age) * e.current.velocity_adjuster;

		//----Velocity
		tfxVec2 current_velocity;
		if (p.velocity_scale) {
			if (e.properties.flags & tfxEmitterPropertyFlags_relative_position) {
				current_velocity = ((p.base.velocity * p.velocity_scale) + p.motion_randomness_speed) * velocity_normal * UPDATE_TIME;
				current_velocity.y += p.weight_acceleration * UPDATE_TIME;
			}
			else {
				current_velocity = ((p.base.velocity * p.velocity_scale) + p.motion_randomness_speed) * velocity_normal * UPDATE_TIME * e.world.scale;
				current_velocity.y += p.weight_acceleration * UPDATE_TIME * e.world.scale.y;
			}
		}

		//Direction
		if (e.properties.emission_type == tfxLine && e.properties.flags & tfxEmitterPropertyFlags_edge_traversal) {
			tfxVec2 offset = velocity_normal * e.current.emitter_size.y;
			float length = std::fabsf(std::sqrtf(current_velocity.Squared()));
			p.distance_travelled += length;

			float emitter_length = e.current.emitter_size.y;
			if (p.distance_travelled > emitter_length) {
				if (e.properties.end_behaviour == LineTraversalEndBehaviour::tfxLoop) {
					p.local.position = p.local.position - offset;
					p.distance_travelled -= emitter_length;
					TransformParticle(p, e);
					p.captured.position = p.world.position;
				}
				else if (e.properties.end_behaviour == LineTraversalEndBehaviour::tfxKill) {
					return false;
				}
				else if (e.properties.end_behaviour == LineTraversalEndBehaviour::tfxLetFree) {
				}
			}
		}
		else {

		}

		//----Color changes
		p.color.a = unsigned char(255.f * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].opacity, p.age, p.max_age) * e.parent->current.color.a);
		p.intensity = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].intensity, p.age, p.max_age);
		if (!(e.properties.flags & tfxEmitterPropertyFlags_random_color)) {
			p.color.r = unsigned char(255.f * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].red, p.age, p.max_age));
			p.color.g = unsigned char(255.f * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].green, p.age, p.max_age));
			p.color.b = unsigned char(255.f * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].blue, p.age, p.max_age));
		}

		p.color = tfxRGBA8(p.color.r, p.color.g, p.color.b, p.color.a);

		//----Size Changes
		p.local.scale.x = p.base.size.x * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].width, p.age, p.max_age);
		//Just here to test:
		//float test1 = p.base.size.x * e.library->LookupPreciseOvertimeNodeList(tfxOvertime_width, e.lookup_node_index, p.age, p.max_age);
		//p.local.scale.x = p.base.size.x * e.library->LookupFastOvertimeValueList(tfxOvertime_width, e.lookup_value_index, p.age, p.max_age);
		if (p.local.scale.x < 0.f)
			p.local.scale.x = p.local.scale.x;

		//----Stretch Changes
		float stretch = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].stretch, p.age, p.max_age);
		if (stretch) {
			float velocity = std::fabsf(p.velocity_scale * p.base.velocity + p.motion_randomness_speed + p.weight_acceleration);
			if (e.properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size)
				p.local.scale.y = (lookup_overtime_callback(e.library->overtime_graphs[e.overtime].width, p.age, p.max_age) *
				(p.base.height + (velocity * stretch * e.current.stretch))) / e.properties.image->image_size.y;
			else
				p.local.scale.y = (lookup_overtime_callback(e.library->overtime_graphs[e.overtime].height, p.age, p.max_age) *
				(p.base.height + (velocity * stretch * e.current.stretch))) / e.properties.image->image_size.y;
			if (p.local.scale.y < p.local.scale.x)
				p.local.scale.y = p.local.scale.x;
		}
		else {
			if (e.properties.flags & tfxEmitterPropertyFlags_lifetime_uniform_size) {
				p.local.scale.y = p.base.size.y * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].width, p.age, p.max_age);
			}
			else {
				p.local.scale.y = p.base.size.y * lookup_overtime_callback(e.library->overtime_graphs[e.overtime].height, p.age, p.max_age);
			}
		}

		//----Spin and angle Changes
		float spin = 0;
		if (e.properties.angle_setting != AngleSetting::tfxAlign && !(e.properties.flags & tfxEmitterPropertyFlags_relative_angle)) {
			float test = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].spin, p.age, p.max_age);
			spin = p.base.spin * test;
		}

		//----Image Changes
		p.image_frame_rate = lookup_overtime_callback(e.library->overtime_graphs[e.overtime].frame_rate, p.age, p.max_age);

		//---------------
		//Now that the latest changes are applied, affect the particle state
		//---------------

		//Before applying the behaviour to the particle position, scale and rotation, you have the chance to override them here
		if (e.particle_update_callback)
			e.particle_update_callback(p);

		//----Rotation
		if (e.properties.angle_setting != AngleSetting::tfxAlign && !(e.properties.flags & tfxEmitterPropertyFlags_relative_angle)) {
			if (spin)
				p.local.rotation += spin * UPDATE_TIME;
		}
		else if (e.properties.angle_setting == AngleSetting::tfxAlign) {
			p.local.rotation = GetVectorAngle(velocity_normal.x, velocity_normal.y) + e.properties.angle_offset;
		}

		//----Position
		p.local.position += current_velocity;

		//----Image animation
		if (e.properties.flags & tfxEmitterPropertyFlags_animate) {
			if (e.properties.flags & tfxEmitterPropertyFlags_reverse_animation)
				p.image_frame -= p.image_frame_rate * UPDATE_TIME;
			else
				p.image_frame += p.image_frame_rate * UPDATE_TIME;

			if (p.image_frame >= e.properties.end_frame + 1) {
				if (e.properties.flags & tfxEmitterPropertyFlags_play_once)
					p.image_frame = e.properties.end_frame;
				else
					p.image_frame -= e.properties.end_frame + 1;
			}
			else if (p.image_frame < 0) {
				if (e.properties.flags & tfxEmitterPropertyFlags_play_once)
					p.image_frame = 0;
				else
					p.image_frame += e.properties.end_frame;
			}
		}
		return true;
	}

	EffectEmitter CreateEffector(float x, float y) {
		EffectEmitter effector;
		effector.local.position = tfxVec2(x, y);

		return effector;
	}

	void EffectEmitter::TransformEffector(EffectEmitter &parent, bool relative_position, bool relative_angle) {
		float s = sin(local.rotation);
		float c = cos(local.rotation);
		matrix.Set(c, s, -s, c);

		if (relative_position) {
			world.scale = local.scale;

			world.rotation = parent.world.rotation + local.rotation;

			matrix = matrix.Transform(parent.matrix);
			tfxVec2 rotatevec = parent.matrix.TransformVector(tfxVec2(local.position.x, local.position.y));

			world.position = parent.world.position + rotatevec * parent.world.scale;

		}
		else {
			world = local;
		}
	}

	void EffectEmitter::TransformEffector(Particle &parent, bool relative_position, bool relative_angle) {
		float s = sin(local.rotation);
		float c = cos(local.rotation);

		matrix.Set(c, s, -s, c);

		if (relative_position) {
			world.scale = local.scale;

			world.rotation = parent.world.rotation + local.rotation;

			matrix = matrix.Transform(parent.matrix);
			tfxVec2 rotatevec = parent.matrix.TransformVector(tfxVec2(local.position.x, local.position.y));

			world.position = parent.world.position + rotatevec * parent.world.scale;

		}
		else {
			world = local;
		}

	}

	void TransformParticle(Particle &p, EffectEmitter &e) {
		//The Particle matrix is only needed for sub effect transformations
		float s = sin(p.local.rotation);
		float c = cos(p.local.rotation);
		p.matrix.Set(c, s, -s, c);
		bool line = (e.properties.flags & tfxEmitterPropertyFlags_edge_traversal && e.properties.emission_type == tfxLine);

		if (e.properties.flags & tfxEmitterPropertyFlags_relative_position || line) {
			p.world.scale = p.local.scale;

			if (e.properties.flags & tfxEmitterPropertyFlags_relative_angle || line)
				p.world.rotation = e.world.rotation + p.local.rotation;
			else
				p.world.rotation = p.local.rotation;

			p.matrix = p.matrix.Transform(e.matrix);
			tfxVec2 rotatevec = e.matrix.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));

			p.world.position = e.world.position + rotatevec * e.world.scale;

		}
		else {
			p.world.position = p.local.position;
			p.world.scale = p.local.scale;
			if (e.properties.flags & tfxEmitterPropertyFlags_relative_angle)
				p.world.rotation = e.world.rotation + p.local.rotation;
			else
				p.world.rotation = p.local.rotation;
		}

	}

	void TransformParticlePrevious(Particle &p, EffectEmitter &e) {
		//The point of this function is so that newly spawned particles can be spawned at the correct tween point
		//p.matrix.Set(std::cosf(p.local.rotation), std::sinf(p.local.rotation), -std::sinf(p.local.rotation), std::cosf(p.local.rotation));
		bool line = (e.properties.flags & tfxEmitterPropertyFlags_edge_traversal && e.properties.emission_type == tfxLine);

		if (e.properties.flags & tfxEmitterPropertyFlags_relative_position || line) {
			p.world.scale = p.local.scale * e.captured.scale;

			if (e.properties.flags & tfxEmitterPropertyFlags_relative_angle || line || e.properties.angle_setting == AngleSetting::tfxAlign)
				p.world.rotation = e.captured.rotation + p.local.rotation;
			else
				p.world.rotation = p.local.rotation;

			Matrix2 mat;
			float s = sin(e.captured.rotation);
			float c = cos(e.captured.rotation);
			mat.Set(c, s, -s, c);
			//mat = mat.Transform(e.matrix);

			//p.matrix = p.matrix.Transform(mat);

			tfxVec2 rotatevec = mat.TransformVector(tfxVec2(p.local.position.x, p.local.position.y));

			p.world.position = e.captured.position + rotatevec * e.captured.scale;

		}
		else {
			p.world.position = p.local.position;
			p.world.scale = p.local.scale;
			if (e.properties.flags & tfxEmitterPropertyFlags_relative_angle || e.properties.angle_setting == AngleSetting::tfxAlign)
				p.world.rotation = e.world.rotation + p.local.rotation;
			else
				p.world.rotation = p.local.rotation;
		}

	}

	void EffectEmitter::ResetGlobalGraphs(bool add_node) {
		library->global_graphs[global].life.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].life.type = tfxGlobal_life;
		library->global_graphs[global].amount.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].amount.type = tfxGlobal_amount;
		library->global_graphs[global].velocity.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].velocity.type = tfxGlobal_velocity;
		library->global_graphs[global].width.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].width.type = tfxGlobal_width;
		library->global_graphs[global].height.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].height.type = tfxGlobal_height;
		library->global_graphs[global].weight.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].weight.type = tfxGlobal_weight;
		library->global_graphs[global].spin.Reset(1.f, tfxGlobalPercentPresetSigned, add_node); library->global_graphs[global].spin.type = tfxGlobal_spin;
		library->global_graphs[global].effect_angle.Reset(0.f, tfxAnglePreset, add_node); library->global_graphs[global].effect_angle.type = tfxGlobal_effect_angle;
		library->global_graphs[global].stretch.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].stretch.type = tfxGlobal_stretch;
		library->global_graphs[global].overal_scale.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].overal_scale.type = tfxGlobal_overal_scale;
		library->global_graphs[global].opacity.Reset(1.f, tfxGlobalOpacityPreset, add_node); library->global_graphs[global].opacity.type = tfxGlobal_opacity;
		library->global_graphs[global].frame_rate.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].frame_rate.type = tfxGlobal_frame_rate;
		library->global_graphs[global].splatter.Reset(1.f, tfxGlobalPercentPreset, add_node); library->global_graphs[global].splatter.type = tfxGlobal_splatter;
		library->CompileGlobalGraph(global);
	}

	void EffectEmitter::ResetBaseGraphs(bool add_node) {
		library->base_graphs[base].life.Reset(1000.f, tfxLifePreset, add_node); library->base_graphs[base].life.type = tfxBase_life;
		library->base_graphs[base].amount.Reset(1.f, tfxAmountPreset, add_node); library->base_graphs[base].amount.type = tfxBase_amount;
		library->base_graphs[base].velocity.Reset(0.f, tfxVelocityPreset, add_node); library->base_graphs[base].velocity.type = tfxBase_velocity;
		library->base_graphs[base].width.Reset(128.f, tfxDimensionsPreset, add_node); library->base_graphs[base].width.type = tfxBase_width;
		library->base_graphs[base].height.Reset(128.f, tfxDimensionsPreset, add_node); library->base_graphs[base].height.type = tfxBase_height;
		library->base_graphs[base].weight.Reset(0.f, tfxWeightPreset, add_node); library->base_graphs[base].weight.type = tfxBase_weight;
		library->base_graphs[base].spin.Reset(0.f, tfxSpinPreset, add_node); library->base_graphs[base].spin.type = tfxBase_spin;
		library->CompileBaseGraph(base);
	}

	void EffectEmitter::ResetPropertyGraphs(bool add_node) {
		library->property_graphs[property].emission_angle.Reset(0.f, tfxAnglePreset, add_node); library->property_graphs[property].emission_angle.type = tfxProperty_emission_angle;
		library->property_graphs[property].emission_range.Reset(0.f, tfxEmissionRangePreset, add_node); library->property_graphs[property].emission_range.type = tfxProperty_emission_range;
		library->property_graphs[property].emitter_angle.Reset(0.f, tfxAnglePreset, add_node); library->property_graphs[property].emitter_angle.type = tfxProperty_emitter_angle;
		library->property_graphs[property].splatter.Reset(0.f, tfxDimensionsPreset, add_node); library->property_graphs[property].splatter.type = tfxProperty_splatter;
		library->property_graphs[property].emitter_width.Reset(0.f, tfxDimensionsPreset, add_node); library->property_graphs[property].emitter_width.type = tfxProperty_emitter_width;
		library->property_graphs[property].emitter_height.Reset(0.f, tfxDimensionsPreset, add_node); library->property_graphs[property].emitter_height.type = tfxProperty_emitter_height;
		library->property_graphs[property].arc_size.Reset(tfxRadians(360.f), tfxArcPreset, add_node); library->property_graphs[property].arc_size.type = tfxProperty_arc_size;
		library->property_graphs[property].arc_offset.Reset(0.f, tfxArcPreset, add_node); library->property_graphs[property].arc_offset.type = tfxProperty_arc_offset;
		library->CompilePropertyGraph(property);
	}

	void EffectEmitter::ResetVariationGraphs(bool add_node) {
		library->variation_graphs[variation].life.Reset(0.f, tfxLifePreset, add_node); library->variation_graphs[variation].life.type = tfxVariation_life;
		library->variation_graphs[variation].amount.Reset(0.f, tfxAmountPreset, add_node); library->variation_graphs[variation].amount.type = tfxVariation_amount;
		library->variation_graphs[variation].velocity.Reset(0.f, tfxVelocityPreset, add_node); library->variation_graphs[variation].velocity.type = tfxVariation_velocity;
		library->variation_graphs[variation].width.Reset(0.f, tfxDimensionsPreset, add_node); library->variation_graphs[variation].width.type = tfxVariation_width;
		library->variation_graphs[variation].height.Reset(0.f, tfxDimensionsPreset, add_node); library->variation_graphs[variation].height.type = tfxVariation_height;
		library->variation_graphs[variation].weight.Reset(0.f, tfxWeightVariationPreset, add_node); library->variation_graphs[variation].weight.type = tfxVariation_weight;
		library->variation_graphs[variation].spin.Reset(0.f, tfxSpinVariationPreset, add_node); library->variation_graphs[variation].spin.type = tfxVariation_spin;
		library->variation_graphs[variation].motion_randomness.Reset(0.f, tfxGlobalPercentPreset, add_node); library->variation_graphs[variation].motion_randomness.type = tfxVariation_motion_randomness;
		library->CompileVariationGraph(variation);
	}

	void EffectEmitter::ResetOvertimeGraphs(bool add_node) {
		library->overtime_graphs[overtime].velocity.Reset(1.f, tfxVelocityOvertimePreset, add_node); library->overtime_graphs[overtime].velocity.type = tfxOvertime_velocity;
		library->overtime_graphs[overtime].velocity_adjuster.Reset(1.f, tfxGlobalPercentPreset, add_node); library->overtime_graphs[overtime].velocity_adjuster.type = tfxOvertime_velocity_adjuster;
		library->overtime_graphs[overtime].width.Reset(1.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].width.type = tfxOvertime_width;
		library->overtime_graphs[overtime].height.Reset(1.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].height.type = tfxOvertime_height;
		library->overtime_graphs[overtime].weight.Reset(1.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].weight.type = tfxOvertime_weight;
		library->overtime_graphs[overtime].spin.Reset(0.f, tfxSpinOvertimePreset, add_node); library->overtime_graphs[overtime].spin.type = tfxOvertime_spin;
		library->overtime_graphs[overtime].stretch.Reset(0.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].stretch.type = tfxOvertime_stretch;
		library->overtime_graphs[overtime].red.Reset(1.f, tfxColorPreset, add_node); library->overtime_graphs[overtime].red.type = tfxOvertime_red;
		library->overtime_graphs[overtime].green.Reset(1.f, tfxColorPreset, add_node); library->overtime_graphs[overtime].green.type = tfxOvertime_green;
		library->overtime_graphs[overtime].blue.Reset(1.f, tfxColorPreset, add_node); library->overtime_graphs[overtime].blue.type = tfxOvertime_blue;
		library->overtime_graphs[overtime].opacity.Reset(1.f, tfxOpacityOvertimePreset, add_node); library->overtime_graphs[overtime].opacity.type = tfxOvertime_opacity;
		library->overtime_graphs[overtime].intensity.Reset(1.f, tfxIntensityOvertimePreset, add_node); library->overtime_graphs[overtime].intensity.type = tfxOvertime_intensity;
		library->overtime_graphs[overtime].frame_rate.Reset(30.f, tfxFrameratePreset, add_node); library->overtime_graphs[overtime].frame_rate.type = tfxOvertime_frame_rate;
		library->overtime_graphs[overtime].stretch.Reset(0.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].stretch.type = tfxOvertime_stretch;
		library->overtime_graphs[overtime].motion_randomness.Reset(0.f, tfxPercentOvertime, add_node); library->overtime_graphs[overtime].motion_randomness.type = tfxOvertime_motion_randomness;
		library->overtime_graphs[overtime].direction.Reset(0.f, tfxDirectionOvertimePreset, add_node); library->overtime_graphs[overtime].direction.type = tfxOvertime_direction;
		library->CompileOvertimeGraph(overtime);
	}

	void EffectEmitter::ResetEffectGraphs(bool add_node) {
		ResetGlobalGraphs(add_node);
	}

	void EffectEmitter::ResetEmitterGraphs(bool add_node) {
		ResetBaseGraphs(add_node);
		ResetPropertyGraphs(add_node);
		ResetVariationGraphs(add_node);
		UpdateMaxLife();
		ResetOvertimeGraphs(add_node);
	}

	void EffectEmitter::InitialiseUninitialisedGraphs() {
		if (type == tfxEffect) {
			if (library->global_graphs[global].life.nodes.size() == 0) library->global_graphs[global].life.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].amount.nodes.size() == 0) library->global_graphs[global].amount.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].velocity.nodes.size() == 0) library->global_graphs[global].velocity.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].width.nodes.size() == 0) library->global_graphs[global].width.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].height.nodes.size() == 0) library->global_graphs[global].height.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].weight.nodes.size() == 0) library->global_graphs[global].weight.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].spin.nodes.size() == 0) library->global_graphs[global].spin.Reset(1.f, tfxGlobalPercentPresetSigned);
			if (library->global_graphs[global].effect_angle.nodes.size() == 0) library->global_graphs[global].effect_angle.Reset(0.f, tfxAnglePreset);
			if (library->global_graphs[global].stretch.nodes.size() == 0) library->global_graphs[global].stretch.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].overal_scale.nodes.size() == 0) library->global_graphs[global].overal_scale.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].opacity.nodes.size() == 0) library->global_graphs[global].opacity.Reset(1.f, tfxOpacityOvertimePreset);
			if (library->global_graphs[global].frame_rate.nodes.size() == 0) library->global_graphs[global].frame_rate.Reset(1.f, tfxGlobalPercentPreset);
			if (library->global_graphs[global].splatter.nodes.size() == 0) library->global_graphs[global].splatter.Reset(1.f, tfxGlobalPercentPreset);
		}

		if (type == tfxEmitter) {
			if (library->base_graphs[base].life.nodes.size() == 0) library->base_graphs[base].life.Reset(1000.f, tfxLifePreset);
			if (library->base_graphs[base].amount.nodes.size() == 0) library->base_graphs[base].amount.Reset(1.f, tfxAmountPreset);
			if (library->base_graphs[base].velocity.nodes.size() == 0) library->base_graphs[base].velocity.Reset(0.f, tfxVelocityPreset);
			if (library->base_graphs[base].width.nodes.size() == 0) library->base_graphs[base].width.Reset(128.f, tfxDimensionsPreset);
			if (library->base_graphs[base].height.nodes.size() == 0) library->base_graphs[base].height.Reset(128.f, tfxDimensionsPreset);
			if (library->base_graphs[base].weight.nodes.size() == 0) library->base_graphs[base].weight.Reset(0.f, tfxWeightPreset);
			if (library->base_graphs[base].spin.nodes.size() == 0) library->base_graphs[base].spin.Reset(0.f, tfxSpinPreset);

			if (library->property_graphs[property].emitter_angle.nodes.size() == 0) library->property_graphs[property].emitter_angle.Reset(0.f, tfxAnglePreset);
			if (library->property_graphs[property].emission_angle.nodes.size() == 0) library->property_graphs[property].emission_angle.Reset(0.f, tfxAnglePreset);
			if (library->property_graphs[property].emission_range.nodes.size() == 0) library->property_graphs[property].emission_range.Reset(0.f, tfxEmissionRangePreset);
			if (library->property_graphs[property].splatter.nodes.size() == 0) library->property_graphs[property].splatter.Reset(0.f, tfxDimensionsPreset);
			if (library->property_graphs[property].emitter_width.nodes.size() == 0) library->property_graphs[property].emitter_width.Reset(0.f, tfxDimensionsPreset);
			if (library->property_graphs[property].emitter_height.nodes.size() == 0) library->property_graphs[property].emitter_height.Reset(0.f, tfxDimensionsPreset);
			if (library->property_graphs[property].arc_size.nodes.size() == 0) library->property_graphs[property].arc_size.Reset(tfxRadians(360.f), tfxArcPreset);
			if (library->property_graphs[property].arc_offset.nodes.size() == 0) library->property_graphs[property].arc_offset.Reset(0.f, tfxArcPreset);

			if (library->variation_graphs[variation].life.nodes.size() == 0) library->variation_graphs[variation].life.Reset(0.f, tfxLifePreset);
			if (library->variation_graphs[variation].amount.nodes.size() == 0) library->variation_graphs[variation].amount.Reset(0.f, tfxAmountPreset);
			if (library->variation_graphs[variation].velocity.nodes.size() == 0) library->variation_graphs[variation].velocity.Reset(0.f, tfxVelocityPreset);
			if (library->variation_graphs[variation].weight.nodes.size() == 0) library->variation_graphs[variation].weight.Reset(0.f, tfxVelocityPreset);
			if (library->variation_graphs[variation].width.nodes.size() == 0) library->variation_graphs[variation].width.Reset(0.f, tfxDimensionsPreset);
			if (library->variation_graphs[variation].height.nodes.size() == 0) library->variation_graphs[variation].height.Reset(0.f, tfxDimensionsPreset);
			if (library->variation_graphs[variation].weight.nodes.size() == 0) library->variation_graphs[variation].weight.Reset(0.f, tfxWeightVariationPreset);
			if (library->variation_graphs[variation].spin.nodes.size() == 0) library->variation_graphs[variation].spin.Reset(0.f, tfxSpinVariationPreset);
			if (library->variation_graphs[variation].motion_randomness.nodes.size() == 0) library->variation_graphs[variation].motion_randomness.Reset(0.f, tfxGlobalPercentPreset);

			if (library->overtime_graphs[overtime].velocity.nodes.size() == 0) library->overtime_graphs[overtime].velocity.Reset(1.f, tfxVelocityOvertimePreset);
			if (library->overtime_graphs[overtime].width.nodes.size() == 0) library->overtime_graphs[overtime].width.Reset(1.f, tfxPercentOvertime);
			if (library->overtime_graphs[overtime].height.nodes.size() == 0) library->overtime_graphs[overtime].height.Reset(1.f, tfxPercentOvertime);
			if (library->overtime_graphs[overtime].weight.nodes.size() == 0) library->overtime_graphs[overtime].weight.Reset(1.f, tfxPercentOvertime);
			if (library->overtime_graphs[overtime].spin.nodes.size() == 0) library->overtime_graphs[overtime].spin.Reset(0.f, tfxSpinOvertimePreset);
			if (library->overtime_graphs[overtime].stretch.nodes.size() == 0) library->overtime_graphs[overtime].stretch.Reset(0.f, tfxPercentOvertime);
			if (library->overtime_graphs[overtime].red.nodes.size() == 0) library->overtime_graphs[overtime].red.Reset(1.f, tfxColorPreset);
			if (library->overtime_graphs[overtime].green.nodes.size() == 0) library->overtime_graphs[overtime].green.Reset(1.f, tfxColorPreset);
			if (library->overtime_graphs[overtime].blue.nodes.size() == 0) library->overtime_graphs[overtime].blue.Reset(1.f, tfxColorPreset);
			if (library->overtime_graphs[overtime].opacity.nodes.size() == 0) library->overtime_graphs[overtime].opacity.Reset(1.f, tfxOpacityOvertimePreset);
			if (library->overtime_graphs[overtime].intensity.nodes.size() == 0) library->overtime_graphs[overtime].intensity.Reset(1.f, tfxIntensityOvertimePreset);
			if (library->overtime_graphs[overtime].frame_rate.nodes.size() == 0) library->overtime_graphs[overtime].frame_rate.Reset(30.f, tfxFrameratePreset);
			if (library->overtime_graphs[overtime].motion_randomness.nodes.size() == 0) library->overtime_graphs[overtime].motion_randomness.Reset(0.f, tfxPercentOvertime);
			if (library->overtime_graphs[overtime].velocity_adjuster.nodes.size() == 0) library->overtime_graphs[overtime].velocity_adjuster.Reset(1.f, tfxGlobalPercentPreset);
			if (library->overtime_graphs[overtime].direction.nodes.size() == 0) library->overtime_graphs[overtime].direction.Reset(0.f, tfxDirectionOvertimePreset);
		}
	}

	void EffectEmitter::SetName(const char *n) {
		name = n;
	}

	FormState Tween(float tween, FormState &world, FormState &captured) {
		FormState tweened;
		tweened.position = world.position * tween + captured.position * (1.f - tween);
		tweened.scale = world.scale * tween + captured.scale * (1.f - tween);
		//Not tweening rotation for now, need to figure out when it tweens over 180 degrees.
		//tweened.rotation = world.rotation * tween + captured.rotation * (1.f - tween);
		tweened.rotation = world.rotation;

		return tweened;
	}

	void EffectEmitter::ClearColors() {
		library->overtime_graphs[overtime].red.Clear();
		library->overtime_graphs[overtime].green.Clear();
		library->overtime_graphs[overtime].blue.Clear();
	}

	void EffectEmitter::AddColorOvertime(float frame, tfxRGB color) {
		library->overtime_graphs[overtime].red.AddNode(frame, color.r);
		library->overtime_graphs[overtime].green.AddNode(frame, color.g);
		library->overtime_graphs[overtime].blue.AddNode(frame, color.b);
	}

	void EffectEmitter::ReSeed(uint64_t seed) {
		if (seed == 0) {
			seed = 0xFFFFFFFFFFF;
		}
		random_generation.ReSeed(seed, seed / 2);
	}

	 void EffectEmitter::SetUpdateCallback(void(*callback)(EffectEmitter &effectemitter)) {
		update_callback = callback;
	}

	 void EffectEmitter::SetUserData(void *data) {
		user_data = data;
	}

	void* EffectEmitter::GetUserData() {
		return user_data;
	}

	void EffectEmitter::SetTimeout(unsigned int frames) {
		timeout = frames;
		for (auto &sub : sub_effectors) {
			sub.SetTimeout(frames);
		}
	}

	bool EffectEmitter::HasSingle() {
		for (auto &e : sub_effectors) {
			if (e.properties.flags & (tfxEmitterPropertyFlags_single | tfxEmitterPropertyFlags_one_shot))
				return true;
		}
		return false;
	}

	bool EffectEmitter::RenameSubEffector(EffectEmitter &emitter, const char *new_name) {
		if (!NameExists(emitter, new_name) && strlen(new_name) > 0) {
			emitter.SetName(new_name);
			library->UpdateEffectPaths();
			return true;
		}

		return false;
	}

	bool EffectEmitter::NameExists(EffectEmitter &emitter, const char *name) {
		for (auto &e : sub_effectors) {
			if (&emitter != &e) {
				if (e.name == name) {
					return true;
				}
			}
		}

		return false;
	}

	void EffectEmitter::ReIndex() {
		unsigned int index = 0;
		for (auto &e : sub_effectors) {
			e.library_index = index++;
			e.parent = this;
			e.ReIndex();
		}
	}

	EffectEmitter* EffectEmitter::GetRootEffect() {
		if (!parent)
			return nullptr;
		EffectEmitter *p = parent;
		unsigned int timeout = 0;
		while (p || ++timeout < 100) {
			if (!p->parent)
				return p;
			p = p->parent;
		}
		return nullptr;
	}

	void EffectEmitter::ResetParents() {
		parent = nullptr;
		root_parent = nullptr;
		for (auto &e : sub_effectors) {
			e.ResetParents();
		}
	}

	EffectEmitter* EffectEmitter::MoveUp(EffectEmitter &emitter) {
		if (emitter.library_index > 0) {
			unsigned int new_index = emitter.library_index - 1;
			std::swap(sub_effectors[emitter.library_index], sub_effectors[new_index]);
			ReIndex();
			emitter.library->UpdateEffectPaths();
			return &sub_effectors[new_index];
		}

		return nullptr;
	}

	EffectEmitter* EffectEmitter::MoveDown(EffectEmitter &emitter) {
		if (emitter.library_index < sub_effectors.size() - 1) {
			unsigned int new_index = emitter.library_index + 1;
			std::swap(sub_effectors[emitter.library_index], sub_effectors[new_index]);
			ReIndex();
			emitter.library->UpdateEffectPaths();
			return &sub_effectors[new_index];
		}
		return nullptr;
	}

	void EffectEmitter::DeleteEmitter(EffectEmitter *emitter) {
		tfxvec<EffectEmitter> stack;
		stack.push_back(*emitter);
		while (stack.size()) {
			EffectEmitter &current = stack.pop_back();
			if (current.type == tfxEffect && !current.parent) {
				library->FreeGlobal(current.global);
			}
			else if(current.type == tfxEmitter) {
				library->FreeProperty(current.property);
				library->FreeBase(current.base);
				library->FreeVariation(current.variation);
				library->FreeOvertime(current.overtime);
			}
			for (auto &sub : current.sub_effectors) {
				stack.push_back(sub);
			}
		}
		sub_effectors.erase(emitter);

		ReIndex();
	}

	void EffectEmitter::CleanUp() {
		if (sub_effectors.size()) {
			tfxvec<EffectEmitter> stack;
			stack.push_back(*this);
			while (stack.size()) {
				EffectEmitter current = stack.pop_back();
				if (current.type == tfxEffect) {
					library->FreeGlobal(current.global);
				}
				else {
					library->FreeProperty(current.property);
					library->FreeBase(current.base);
					library->FreeVariation(current.variation);
					library->FreeOvertime(current.overtime);
				}
				for (auto &sub : current.sub_effectors) {
					stack.push_back(sub);
				}
				current.sub_effectors.clear();
			}
		}

		ReIndex();
	}

	void EffectEmitter::Clone(EffectEmitter &clone, EffectEmitter *root_parent, EffectLibrary *destination_library) {
		clone = *this;
		clone.sub_effectors.clear();
		clone.flags |= tfxEmitterStateFlags_enabled;
		clone.user_data = nullptr;
		clone.library = destination_library;

		if (type == tfxEffect) {
			if (root_parent == &clone) {
				clone.global = library->CloneGlobal(global, destination_library);
				clone.library->CompileGlobalGraph(clone.global);
			}
			else {
				clone.global = root_parent->global;
			}
		}
		else {
			clone.property = library->CloneProperty(property, destination_library);
			clone.base = library->CloneBase(base, destination_library);
			clone.variation = library->CloneVariation(variation, destination_library);
			clone.overtime = library->CloneOvertime(overtime, destination_library);
			clone.UpdateMaxLife();
			clone.library->CompilePropertyGraph(clone.property);
			clone.library->CompileBaseGraph(clone.base);
			clone.library->CompileVariationGraph(clone.variation);
			clone.library->CompileOvertimeGraph(clone.overtime);
		}

		for (auto &e : sub_effectors) {
			if (e.type == tfxEmitter) {
				EffectEmitter emitter_copy;
				e.Clone(emitter_copy, root_parent, destination_library);
				emitter_copy.user_data = nullptr;
				clone.AddEmitter(emitter_copy);
			}
			else {
				EffectEmitter effect_copy;
				e.Clone(effect_copy, root_parent, destination_library);
				effect_copy.user_data = nullptr;
				clone.AddEffect(effect_copy);
			}
		}
	}

	void EffectEmitter::EnableAllEmitters() {
		for (auto &e : sub_effectors) {
			e.flags |= tfxEmitterStateFlags_enabled;
			e.EnableAllEmitters();
		}
	}

	void EffectEmitter::EnableEmitter() {
		flags |= tfxEmitterStateFlags_enabled;
	}

	void EffectEmitter::DisableAllEmitters() {
		for (auto &e : sub_effectors) {
			e.flags &= ~tfxEmitterStateFlags_enabled;
			e.DisableAllEmitters();
		}
	}

	void EffectEmitter::DisableAllEmittersExcept(EffectEmitter &emitter) {
		for (auto &e : sub_effectors) {
			if (e.library_index == emitter.library_index)
				e.flags |= tfxEmitterStateFlags_enabled;
			else
				e.flags &= ~tfxEmitterStateFlags_enabled;
		}
	}

	Graph* EffectEmitter::GetGraphByType(GraphType type) {

		if (type < tfxGlobalCount) {
			return &((Graph*)&library->global_graphs[global])[type];
		}
		else if (type >= tfxPropertyStart && type < tfxBaseStart) {
			int ref = type - tfxPropertyStart;
			return &((Graph*)&library->property_graphs[property])[ref];
		}
		else if (type >= tfxBaseStart && type < tfxVariationStart) {
			int ref = type - tfxBaseStart;
			return &((Graph*)&library->base_graphs[base])[ref];
		}
		else if (type >= tfxVariationStart && type < tfxOvertimeStart) {
			int ref = type - tfxVariationStart;
			return &((Graph*)&library->variation_graphs[variation])[ref];
		}
		else if (type >= tfxOvertimeStart) {
			int ref = type - tfxOvertimeStart;
			return &((Graph*)&library->overtime_graphs[overtime])[ref];
		}

		return nullptr;

	}

	unsigned int EffectEmitter::GetGraphIndexByType(GraphType type) {

		if (type < tfxGlobalCount) {
			return global;
		}
		else if (type >= tfxGlobalCount && type < tfxPropertyCount) {
			return property;
		}
		else if (type >= tfxPropertyCount && type < tfxBaseCount) {
			return base;
		}
		else if (type >= tfxBaseCount && type < tfxVariationCount) {
			return variation;
		}
		else if (type >= tfxOvertimeCount) {
			return overtime;
		}

		assert(0);	//Unable to find a graph of that type (shouldn't happend)

		return 0;

	}

	void EffectEmitter::FreeGraphs() {
		if (type == tfxEffect) {
			library->global_graphs[global].life.Free();
			library->global_graphs[global].amount.Free();
			library->global_graphs[global].velocity.Free();
			library->global_graphs[global].width.Free();
			library->global_graphs[global].height.Free();
			library->global_graphs[global].weight.Free();
			library->global_graphs[global].spin.Free();
			library->global_graphs[global].effect_angle.Free();
			library->global_graphs[global].stretch.Free();
			library->global_graphs[global].overal_scale.Free();
			library->global_graphs[global].opacity.Free();
			library->global_graphs[global].frame_rate.Free();
			library->global_graphs[global].splatter.Free();
		}

		if (type == tfxEmitter) {
			library->property_graphs[property].emission_angle.Free();
			library->property_graphs[property].emission_range.Free();
			library->property_graphs[property].emitter_angle.Free();
			library->property_graphs[property].splatter.Free();
			library->property_graphs[property].emitter_width.Free();
			library->property_graphs[property].emitter_height.Free();
			library->property_graphs[property].arc_size.Free();
			library->property_graphs[property].arc_offset.Free();

			library->base_graphs[base].life.Free();
			library->base_graphs[base].amount.Free();
			library->base_graphs[base].velocity.Free();
			library->base_graphs[base].width.Free();
			library->base_graphs[base].height.Free();
			library->base_graphs[base].weight.Free();
			library->base_graphs[base].spin.Free();

			library->variation_graphs[variation].life.Free();
			library->variation_graphs[variation].amount.Free();
			library->variation_graphs[variation].velocity.Free();
			library->variation_graphs[variation].width.Free();
			library->variation_graphs[variation].height.Free();
			library->variation_graphs[variation].weight.Free();
			library->variation_graphs[variation].spin.Free();
			library->variation_graphs[variation].motion_randomness.Free();

			library->overtime_graphs[overtime].velocity.Free();
			library->overtime_graphs[overtime].width.Free();
			library->overtime_graphs[overtime].height.Free();
			library->overtime_graphs[overtime].weight.Free();
			library->overtime_graphs[overtime].spin.Free();
			library->overtime_graphs[overtime].stretch.Free();
			library->overtime_graphs[overtime].red.Free();
			library->overtime_graphs[overtime].green.Free();
			library->overtime_graphs[overtime].blue.Free();
			library->overtime_graphs[overtime].opacity.Free();
			library->overtime_graphs[overtime].intensity.Free();
			library->overtime_graphs[overtime].frame_rate.Free();
			library->overtime_graphs[overtime].motion_randomness.Free();
			library->overtime_graphs[overtime].velocity_adjuster.Free();
			library->overtime_graphs[overtime].direction.Free();
		}
	}

	void EffectEmitter::CompileGraphs() {
		if (type == tfxEffect) {
			for (unsigned int t = (unsigned int)tfxGlobal_life; t != (unsigned int)tfxProperty_emission_angle; ++t) {
				CompileGraph(*GetGraphByType(GraphType(t)));
			}

			for (auto &emitter : sub_effectors) {
				for (unsigned int t = (unsigned int)tfxProperty_emission_angle; t != (unsigned int)tfxOvertime_velocity; ++t) {
					CompileGraph(*GetGraphByType((GraphType)t));
				}

				for (unsigned int t = (unsigned int)tfxOvertime_velocity; t != (unsigned int)tfxGraphMaxIndex; ++t) {
					CompileGraphOvertime(*emitter.GetGraphByType((GraphType)t));
				}
			}
		}
	}

	EffectEmitter& EffectLibrary::operator[] (uint32_t index) {
		return effects[index];
	}

	bool EffectLibrary::RenameEffect(EffectEmitter &effect, const char *new_name) {
		if (!NameExists(effect, new_name) && strlen(new_name) > 0) {
			effect.SetName(new_name);
			UpdateEffectPaths();
			return true;
		}

		return false;
	}

	bool EffectLibrary::NameExists(EffectEmitter& effect, const char *name) {
		for (auto &e : effects) {
			if (effect.library_index != e.library_index) {
				if (e.name == name) {
					return true;
				}
			}
		}

		return false;
	}

	bool EffectLibrary::NameExists2(EffectEmitter& effect, const char *name) {
		for (auto &e : effects) {
			if (e.name == name) {
				return true;
			}
		}
		return false;
	}

	void EffectLibrary::UpdateEffectPaths() {
		effect_paths.Clear();
		for (auto &e : effects) {
			tfxText path = e.name;
			e.path_hash = XXHash64::hash(path.c_str(), path.Length(), 0);
			AddPath(e, path);
		}
	}

	void EffectLibrary::AddPath(EffectEmitter &effectemitter, tfxText path) {
		effect_paths.Insert(path, &effectemitter);
		for (auto &sub : effectemitter.sub_effectors) {
			tfxText sub_path = path;
			sub_path.Appendf("/%s", sub.name.c_str());
			sub.path_hash = XXHash64::hash(sub_path.c_str(), sub_path.Length(), 0);
			AddPath(sub, sub_path);
		}
	}

	EffectEmitter &EffectLibrary::AddEffect(EffectEmitter &effect) {
		effect.library_index = effects.current_size;
		effect.type = tfxEffect;
		effect.uid = ++uid;
		effect.library = this;
		effects.push_back(effect);
		ReIndex();
		UpdateEffectPaths();
		return effects.back();
	}

	EffectEmitter* EffectLibrary::GetEffect(tfxText path) {
		assert(effect_paths.ValidName(path));
		return effect_paths.At(path);
	}

	void EffectLibrary::PrepareEffectTemplate(tfxText path, EffectEmitterTemplate &effect_template) {
		EffectEmitter *effect = GetEffect(path);
		assert(effect);
		assert(effect->type == tfxEffect);
		effect->Clone(effect_template.effect_template, &effect_template.effect_template, this);
		effect_template.AddPath(effect_template.effect_template, effect_template.effect_template.name);
	}

	void EffectLibrary::ReIndex() {
		unsigned int index = 0;
		for (auto &e : effects) {
			e.library_index = index++;
			e.parent = nullptr;
			e.ReIndex();
		}
	}

	void EffectLibrary::UpdateParticleShapeReferences(tfxvec<EffectEmitter> &effects, unsigned int default_index) {
		for (auto &effect : effects) {
			for (auto &emitter : effect.sub_effectors) {
				if(particle_shapes.ValidIntName(emitter.properties.shape_index)) 
					emitter.properties.image = &particle_shapes.AtInt(emitter.properties.shape_index);
				else
					emitter.properties.image = &particle_shapes.AtInt(default_index);
				UpdateParticleShapeReferences(emitter.sub_effectors, default_index);
			}
		}
	}

	EffectEmitter* EffectLibrary::MoveUp(EffectEmitter &effect) {
		if (effect.library_index > 0) {
			unsigned int new_index = effect.library_index - 1;
			std::swap(effects[effect.library_index], effects[new_index]);
			UpdateEffectPaths();
			ReIndex();
			return &effects[new_index];
		}
		return nullptr;
	}

	EffectEmitter* EffectLibrary::MoveDown(EffectEmitter &effect) {
		if (effect.library_index < effects.size() - 1) {
			unsigned int new_index = effect.library_index + 1;
			std::swap(effects[effect.library_index], effects[new_index]);
			UpdateEffectPaths();
			ReIndex();
			return &effects[new_index];
		}
		return nullptr;
	}

	void EffectLibrary::RemoveShape(unsigned int shape_index) {
		particle_shapes.RemoveInt(shape_index);
		for (auto &m : particle_shapes.map) {
			particle_shapes[m.index].shape_index = (unsigned int)m.key;
		}
	}

	void EffectLibrary::DeleteEffect(EffectEmitter *effect) {
		effects[effect->library_index].CleanUp();
		effects.erase(&effects[effect->library_index]);

		UpdateEffectPaths();
		ReIndex();
	}

	unsigned int EffectLibrary::AddGlobal() {
		if (free_global_graphs.size())
			return free_global_graphs.pop_back();
		GlobalAttributes global;
		global_graphs.push_back(global);
		return global_graphs.size() - 1;
	}
	unsigned int EffectLibrary::AddProperty() {
		if (free_property_graphs.size())
			return free_property_graphs.pop_back();
		PropertyAttributes property;
		property_graphs.push_back(property);
		return property_graphs.size() - 1;
	}
	unsigned int EffectLibrary::AddBase() {
		if (free_base_graphs.size())
			return free_base_graphs.pop_back();
		BaseAttributes base;
		base_graphs.push_back(base);
		return base_graphs.size() - 1;
	}
	unsigned int EffectLibrary::AddVariation() {
		if (free_variation_graphs.size())
			return free_variation_graphs.pop_back();
		VariationAttributes variation;
		variation_graphs.push_back(variation);
		return variation_graphs.size() - 1;
	}
	unsigned int EffectLibrary::AddOvertime() {
		if (free_overtime_graphs.size())
			return free_overtime_graphs.pop_back();
		OvertimeAttributes overtime;
		overtime_graphs.push_back(overtime);
		return overtime_graphs.size() - 1;
	}

	void EffectLibrary::FreeGlobal(unsigned int index) {
		assert(index < global_graphs.size());
		free_global_graphs.push_back(index);
	}
	void EffectLibrary::FreeProperty(unsigned int index) {
		assert(index < property_graphs.size());
		free_property_graphs.push_back(index);
	}
	void EffectLibrary::FreeBase(unsigned int index) {
		assert(index < base_graphs.size());
		free_base_graphs.push_back(index);
	}
	void EffectLibrary::FreeVariation(unsigned int index) {
		assert(index < variation_graphs.size());
		free_variation_graphs.push_back(index);
	}
	void EffectLibrary::FreeOvertime(unsigned int index) {
		assert(index < overtime_graphs.size());
		free_overtime_graphs.push_back(index);
	}

	unsigned int EffectLibrary::CloneGlobal(unsigned int source_index, EffectLibrary *destination_library) {
		unsigned int index = destination_library->AddGlobal();
		destination_library->global_graphs[index] = global_graphs[source_index];
		return index;
	}

	unsigned int EffectLibrary::CloneProperty(unsigned int source_index, EffectLibrary *destination_library) {
		unsigned int index = destination_library->AddProperty();
		destination_library->property_graphs[index] = property_graphs[source_index];
		return index;
	}

	unsigned int EffectLibrary::CloneBase(unsigned int source_index, EffectLibrary *destination_library) {
		unsigned int index = destination_library->AddBase();
		destination_library->base_graphs[index] = base_graphs[source_index];
		return index;
	}

	unsigned int EffectLibrary::CloneVariation(unsigned int source_index, EffectLibrary *destination_library) {
		unsigned int index = destination_library->AddVariation();
		destination_library->variation_graphs[index] = variation_graphs[source_index];
		return index;
	}

	unsigned int EffectLibrary::CloneOvertime(unsigned int source_index, EffectLibrary *destination_library) {
		unsigned int index = destination_library->AddOvertime();
		destination_library->overtime_graphs[index] = overtime_graphs[source_index];
		return index;
	}

	void EffectLibrary::AddEmitterGraphs(EffectEmitter& emitter) {
		emitter.property = AddProperty();
		emitter.base = AddBase();
		emitter.variation = AddVariation();
		emitter.overtime = AddOvertime();
	}

	void EffectLibrary::AddEffectGraphs(EffectEmitter& effect) {
		EffectEmitter *root_effect = effect.GetRootEffect();
		if (!root_effect)
			effect.global = AddGlobal();
		else
			effect.global = root_effect->global;
	}

	unsigned int EffectLibrary::AddAnimationSettings(EffectEmitter& effect) {
		assert(effect.type == tfxEffect);
		AnimationSettings a;
		a.frames = 32;
		a.current_frame = 1;
		a.frame_offset = 0;
		a.position = tfxVec2(0.f, 0.f);
		a.frame_size = tfxVec2(256.f, 256.f);
		a.loop = false;
		a.seamless = false;
		a.seed = 0;
		a.zoom = 1.f;
		a.scale = 1.f;
		a.needs_recording = true;
		a.needs_exporting = 0;
		a.color_option = ExportColorOptions::tfxFullColor;
		a.export_option = ExportOptions::tfxSpriteSheet;
		a.export_with_transparency = true;
		animation_settings.push_back(a);
		effect.animation_settings = animation_settings.size() - 1;
		return effect.animation_settings;
	}

	void EffectLibrary::Clear() {
		for (auto &e : effects) {
			e.FreeGraphs();
		}
		effects.free_all();
		particle_shapes.Clear();
		global_graphs.free_all();
		property_graphs.free_all();
		base_graphs.free_all();
		variation_graphs.free_all();
		overtime_graphs.free_all();

		free_global_graphs.free_all();
		free_property_graphs.free_all();
		free_base_graphs.free_all();
		free_variation_graphs.free_all();
		free_overtime_graphs.free_all();
		uid = 0;
	}

	void EffectLibrary::UpdateAllNodes() {
		unsigned int running_node_index = 0;
		unsigned int running_value_index = 0;
		tfxvec<EffectEmitter*> stack;
		all_nodes.clear();
		node_lookup_indexes.clear();
		compiled_lookup_values.clear();
		compiled_lookup_indexes.clear();
		for (auto &effect : effects) {
			stack.push_back(&effect);
			while(!stack.empty()) {
				EffectEmitter *current = stack.pop_back();
				EffectLookUpData lookup_data;
				EffectLookUpData value_lookup_data;
				memset(&lookup_data, 0, sizeof(EffectLookUpData));
				memset(&value_lookup_data, 0, sizeof(EffectLookUpData));
				if (current->type == tfxEffect) {
					for (int i = 0; i != tfxGlobalCount; ++i) {

						Graph &graph = ((Graph*)(&global_graphs[current->global]))[i];
						GraphLookupIndex &index = ((GraphLookupIndex*)&lookup_data)[i];
						index.start_index = running_node_index;
						index.length = graph.nodes.size();
						index.max_life = graph.lookup.life;
						for (auto &node : graph.nodes) {
							all_nodes.push_back(node);
							running_node_index++;
						}

						GraphLookupIndex &value_index = ((GraphLookupIndex*)&value_lookup_data)[i];
						value_index.start_index = running_value_index;
						value_index.length = graph.lookup.values.size();
						value_index.max_life = graph.lookup.life;
						for (auto value : graph.lookup.values) {
							compiled_lookup_values.push_back(value);
							running_value_index++;
						}

					}
				}
				else if (current->type == tfxEmitter) {

					int offset = tfxGlobalCount;

					for (int i = 0; i != tfxPropertyCount; ++i) {
						Graph &graph = ((Graph*)(&property_graphs[current->property]))[i];
						GraphLookupIndex &index = ((GraphLookupIndex*)&lookup_data)[i + offset];
						index.start_index = running_node_index;
						index.length = graph.nodes.size();
						index.max_life = graph.lookup.life;
						for (auto &node : graph.nodes) {
							all_nodes.push_back(node);
							running_node_index++;
						}

						GraphLookupIndex &value_index = ((GraphLookupIndex*)&value_lookup_data)[i + offset];
						value_index.start_index = running_value_index;
						value_index.length = graph.lookup.values.size();
						value_index.max_life = graph.lookup.life;
						for (auto value : graph.lookup.values) {
							compiled_lookup_values.push_back(value);
							running_value_index++;
						}
					}

					offset += tfxPropertyCount;

					for (int i = 0; i != tfxBaseCount; ++i) {
						Graph &graph = ((Graph*)(&base_graphs[current->base]))[i];
						GraphLookupIndex &index = ((GraphLookupIndex*)&lookup_data)[i + offset];
						index.start_index = running_node_index;
						index.length = graph.nodes.size();
						index.max_life = graph.lookup.life;
						for (auto &node : graph.nodes) {
							all_nodes.push_back(node);
							running_node_index++;
						}

						GraphLookupIndex &value_index = ((GraphLookupIndex*)&value_lookup_data)[i + offset];
						value_index.start_index = running_value_index;
						value_index.length = graph.lookup.values.size();
						value_index.max_life = graph.lookup.life;
						for (auto value : graph.lookup.values) {
							compiled_lookup_values.push_back(value);
							running_value_index++;
						}
					}

					offset += tfxBaseCount;

					for (int i = 0; i != tfxVariationCount; ++i) {
						Graph &graph = ((Graph*)(&variation_graphs[current->variation]))[i];
						GraphLookupIndex &index = ((GraphLookupIndex*)&lookup_data)[i + offset];
						index.start_index = running_node_index;
						index.length = graph.nodes.size();
						index.max_life = graph.lookup.life;
						for (auto &node : graph.nodes) {
							all_nodes.push_back(node);
							running_node_index++;
						}

						GraphLookupIndex &value_index = ((GraphLookupIndex*)&value_lookup_data)[i + offset];
						value_index.start_index = running_value_index;
						value_index.length = graph.lookup.values.size();
						value_index.max_life = graph.lookup.life;
						for (auto value : graph.lookup.values) {
							compiled_lookup_values.push_back(value);
							running_value_index++;
						}
					}

					offset += tfxVariationCount;

					for (int i = 0; i != tfxOvertimeCount; ++i) {
						Graph &graph = ((Graph*)(&overtime_graphs[current->overtime]))[i];
						GraphLookupIndex &index = ((GraphLookupIndex*)&lookup_data)[i + offset];
						index.start_index = running_node_index;
						index.length = graph.nodes.size();
						index.max_life = graph.lookup.life;
						for (auto &node : graph.nodes) {
							all_nodes.push_back(node);
							running_node_index++;
						}

						GraphLookupIndex &value_index = ((GraphLookupIndex*)&value_lookup_data)[i + offset];
						value_index.start_index = running_value_index;
						value_index.length = graph.lookup.values.size();
						value_index.max_life = graph.lookup.life;
						for (auto value : graph.lookup.values) {
							compiled_lookup_values.push_back(value);
							running_value_index++;
						}
					}

				}

				node_lookup_indexes.push_back(lookup_data);
				compiled_lookup_indexes.push_back(value_lookup_data);
				current->lookup_node_index = node_lookup_indexes.size() - 1;
				current->lookup_value_index = compiled_lookup_indexes.size() - 1;
				for (auto &sub : current->sub_effectors) {
					stack.push_back(&sub);
				}
			}
		}
	}

	void EffectLibrary::CompileAllGraphs() {
		for (auto &g : global_graphs) {
			CompileGraph(g.amount);
			CompileGraph(g.effect_angle);
			CompileGraph(g.frame_rate);
			CompileGraph(g.height);
			CompileGraph(g.width);
			CompileGraph(g.life);
			CompileGraph(g.opacity);
			CompileGraph(g.overal_scale);
			CompileGraph(g.spin);
			CompileGraph(g.splatter);
			CompileGraph(g.stretch);
			CompileGraph(g.velocity);
			CompileGraph(g.weight);
		}
		for (auto &g : property_graphs) {
			CompileGraph(g.arc_offset);
			CompileGraph(g.arc_size);
			CompileGraph(g.emission_angle);
			CompileGraph(g.emission_range);
			CompileGraph(g.emitter_angle);
			CompileGraph(g.emitter_width);
			CompileGraph(g.emitter_height);
			CompileGraph(g.splatter);
		}
		for (auto &g : base_graphs) {
			CompileGraph(g.amount);
			CompileGraph(g.width);
			CompileGraph(g.height);
			CompileGraph(g.life);
			CompileGraph(g.spin);
			CompileGraph(g.velocity);
			CompileGraph(g.weight);
		}
		for (auto &g : variation_graphs) {
			CompileGraph(g.amount);
			CompileGraph(g.width);
			CompileGraph(g.height);
			CompileGraph(g.life);
			CompileGraph(g.motion_randomness);
			CompileGraph(g.spin);
			CompileGraph(g.velocity);
			CompileGraph(g.weight);
		}
		for (auto &g : overtime_graphs) {
			CompileGraphOvertime(g.red);
			CompileGraphOvertime(g.green);
			CompileGraphOvertime(g.blue);
			CompileGraphOvertime(g.opacity);
			CompileGraphOvertime(g.intensity);
			CompileGraphOvertime(g.frame_rate);
			CompileGraphOvertime(g.width);
			CompileGraphOvertime(g.height);
			CompileGraphOvertime(g.motion_randomness);
			CompileGraphOvertime(g.spin);
			CompileGraphOvertime(g.stretch);
			CompileGraphOvertime(g.velocity);
			CompileGraph(g.velocity_adjuster);
			CompileGraphOvertime(g.weight);
			CompileGraphOvertime(g.direction);
		}
	}

	void EffectLibrary::CompileGlobalGraph(unsigned int index) {
		GlobalAttributes &g = global_graphs[index];
		CompileGraph(g.amount);
		CompileGraph(g.effect_angle);
		CompileGraph(g.frame_rate);
		CompileGraph(g.height);
		CompileGraph(g.width);
		CompileGraph(g.life);
		CompileGraph(g.opacity);
		CompileGraph(g.overal_scale);
		CompileGraph(g.spin);
		CompileGraph(g.splatter);
		CompileGraph(g.stretch);
		CompileGraph(g.velocity);
		CompileGraph(g.weight);
	}
	void EffectLibrary::CompilePropertyGraph(unsigned int index) {
		PropertyAttributes &g = property_graphs[index];
		CompileGraph(g.arc_offset);
		CompileGraph(g.arc_size);
		CompileGraph(g.emission_angle);
		CompileGraph(g.emission_range);
		CompileGraph(g.emitter_angle);
		CompileGraph(g.emitter_width);
		CompileGraph(g.emitter_height);
		CompileGraph(g.splatter);
	}
	void EffectLibrary::CompileBaseGraph(unsigned int index) {
		BaseAttributes &g = base_graphs[index];
		CompileGraph(g.amount);
		CompileGraph(g.width);
		CompileGraph(g.height);
		CompileGraph(g.life);
		CompileGraph(g.spin);
		CompileGraph(g.velocity);
		CompileGraph(g.weight);
	}
	void EffectLibrary::CompileVariationGraph(unsigned int index) {
		VariationAttributes &g = variation_graphs[index];
		CompileGraph(g.amount);
		CompileGraph(g.width);
		CompileGraph(g.height);
		CompileGraph(g.life);
		CompileGraph(g.motion_randomness);
		CompileGraph(g.spin);
		CompileGraph(g.velocity);
		CompileGraph(g.weight);
	}
	void EffectLibrary::CompileOvertimeGraph(unsigned int index) {
		OvertimeAttributes &g = overtime_graphs[index];
		CompileGraphOvertime(g.red);
		CompileGraphOvertime(g.green);
		CompileGraphOvertime(g.blue);
		CompileGraphOvertime(g.opacity);
		CompileGraphOvertime(g.intensity);
		CompileGraphOvertime(g.frame_rate);
		CompileGraphOvertime(g.width);
		CompileGraphOvertime(g.height);
		CompileGraphOvertime(g.motion_randomness);
		CompileGraphOvertime(g.spin);
		CompileGraphOvertime(g.stretch);
		CompileGraphOvertime(g.velocity);
		CompileGraph(g.velocity_adjuster);
		CompileGraphOvertime(g.weight);
		CompileGraphOvertime(g.direction);
	}
	void EffectLibrary::CompileColorGraphs(unsigned int index) {
		OvertimeAttributes &g = overtime_graphs[index];
		CompileGraphOvertime(g.red);
		CompileGraphOvertime(g.green);
		CompileGraphOvertime(g.blue);
	}

	void EffectLibrary::SetMinMaxData() {
		graph_min_max.clear();
		graph_min_max.create_pool(tfxGraphMaxIndex);

		graph_min_max[tfxGlobal_life] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_amount] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_velocity] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_width] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_height] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_weight] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_spin] = GetMinMaxGraphValues(tfxGlobalPercentPresetSigned);
		graph_min_max[tfxGlobal_effect_angle] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxGlobal_stretch] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_overal_scale] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_opacity] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
		graph_min_max[tfxGlobal_frame_rate] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxGlobal_splatter] = GetMinMaxGraphValues(tfxGlobalPercentPreset);

		graph_min_max[tfxProperty_emitter_angle] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxProperty_emission_angle] = GetMinMaxGraphValues(tfxAnglePreset);
		graph_min_max[tfxProperty_emission_range] = GetMinMaxGraphValues(tfxEmissionRangePreset);
		graph_min_max[tfxProperty_splatter] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxProperty_emitter_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxProperty_emitter_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxProperty_arc_size] = GetMinMaxGraphValues(tfxArcPreset);
		graph_min_max[tfxProperty_arc_offset] = GetMinMaxGraphValues(tfxArcPreset);

		graph_min_max[tfxBase_life] = GetMinMaxGraphValues(tfxLifePreset);
		graph_min_max[tfxBase_amount] = GetMinMaxGraphValues(tfxAmountPreset);
		graph_min_max[tfxBase_velocity] = GetMinMaxGraphValues(tfxVelocityPreset);
		graph_min_max[tfxBase_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxBase_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxBase_weight] = GetMinMaxGraphValues(tfxWeightPreset);
		graph_min_max[tfxBase_spin] = GetMinMaxGraphValues(tfxSpinPreset);

		graph_min_max[tfxVariation_life] = GetMinMaxGraphValues(tfxLifePreset);
		graph_min_max[tfxVariation_amount] = GetMinMaxGraphValues(tfxAmountPreset);
		graph_min_max[tfxVariation_velocity] = GetMinMaxGraphValues(tfxVelocityPreset);
		graph_min_max[tfxVariation_width] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxVariation_height] = GetMinMaxGraphValues(tfxDimensionsPreset);
		graph_min_max[tfxVariation_weight] = GetMinMaxGraphValues(tfxWeightVariationPreset);
		graph_min_max[tfxVariation_spin] = GetMinMaxGraphValues(tfxSpinVariationPreset);
		graph_min_max[tfxVariation_motion_randomness] = GetMinMaxGraphValues(tfxGlobalPercentPreset);

		graph_min_max[tfxOvertime_velocity] = GetMinMaxGraphValues(tfxVelocityOvertimePreset);
		graph_min_max[tfxOvertime_width] = GetMinMaxGraphValues(tfxPercentOvertime);
		graph_min_max[tfxOvertime_height] = GetMinMaxGraphValues(tfxPercentOvertime);
		graph_min_max[tfxOvertime_weight] = GetMinMaxGraphValues(tfxPercentOvertime);
		graph_min_max[tfxOvertime_spin] = GetMinMaxGraphValues(tfxSpinOvertimePreset);
		graph_min_max[tfxOvertime_stretch] = GetMinMaxGraphValues(tfxPercentOvertime);
		graph_min_max[tfxOvertime_red] = GetMinMaxGraphValues(tfxColorPreset);
		graph_min_max[tfxOvertime_green] = GetMinMaxGraphValues(tfxColorPreset);
		graph_min_max[tfxOvertime_blue] = GetMinMaxGraphValues(tfxColorPreset);
		graph_min_max[tfxOvertime_opacity] = GetMinMaxGraphValues(tfxOpacityOvertimePreset);
		graph_min_max[tfxOvertime_intensity] = GetMinMaxGraphValues(tfxIntensityOvertimePreset);
		graph_min_max[tfxOvertime_frame_rate] = GetMinMaxGraphValues(tfxFrameratePreset);
		graph_min_max[tfxOvertime_motion_randomness] = GetMinMaxGraphValues(tfxPercentOvertime);
		graph_min_max[tfxOvertime_velocity_adjuster] = GetMinMaxGraphValues(tfxGlobalPercentPreset);
		graph_min_max[tfxOvertime_direction] = GetMinMaxGraphValues(tfxDirectionOvertimePreset);

	}

	float EffectLibrary::LookupPreciseOvertimeNodeList(GraphType graph_type, int lookup_node_index, float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		GraphLookupIndex &lookup_data = ((GraphLookupIndex*)&node_lookup_indexes[lookup_node_index])[graph_type];
		float min_y = graph_min_max[graph_type].y;
		float max_y = graph_min_max[graph_type].w;
		for (int i = lookup_data.start_index; i != lookup_data.start_index + lookup_data.length; ++i) {
			AttributeNode &a = all_nodes[i];
			float frame = a.frame * life;
			if (age < frame) {
				p = (age - lastf) / (frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, min_y, max_y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = frame - 1;
			lastec = &a;
		}
		return lastv;
	}

	float EffectLibrary::LookupPreciseNodeList(GraphType graph_type, int lookup_node_index, float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		GraphLookupIndex &lookup_data = ((GraphLookupIndex*)&node_lookup_indexes[lookup_node_index])[graph_type];
		float min_y = graph_min_max[graph_type].y;
		float max_y = graph_min_max[graph_type].w;
		for (int i = lookup_data.start_index; i != lookup_data.start_index + lookup_data.length; ++i) {
			AttributeNode &a = all_nodes[i];
			if (age < a.frame) {
				p = (age - lastf) / (a.frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, min_y, max_y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = a.frame - 1;
			lastec = &a;
		}
		return lastv;
	}

	float EffectLibrary::LookupFastValueList(GraphType graph_type, int lookup_node_index, float frame) {
		GraphLookupIndex &lookup_data = ((GraphLookupIndex*)&compiled_lookup_indexes[lookup_node_index])[graph_type];
		frame += lookup_data.start_index;
		if ((unsigned int)frame < lookup_data.start_index + lookup_data.length - 1)
			return compiled_lookup_values[(unsigned int)frame];
		return compiled_lookup_values[lookup_data.start_index + lookup_data.length - 1];
	}

	float EffectLibrary::LookupFastOvertimeValueList(GraphType graph_type, int lookup_value_index, float age, float lifetime) {
		GraphLookupIndex &lookup_data = ((GraphLookupIndex*)&compiled_lookup_indexes[lookup_value_index])[graph_type];
		float frame = (float)lookup_data.start_index;
		if (lifetime)
			frame += (age / lifetime * lookup_data.max_life) / tfxLOOKUP_FREQUENCY_OVERTIME;
		if (frame < lookup_data.start_index + lookup_data.length - 1)
			return compiled_lookup_values[(unsigned int)frame];
		return compiled_lookup_values[lookup_data.start_index + lookup_data.length - 1];
	}

	unsigned int EffectLibrary::CountOfGraphsInUse() {
		return global_graphs.size() + property_graphs.size() + base_graphs.size() + variation_graphs.size() + overtime_graphs.size() - CountOfFreeGraphs();
	}

	unsigned int EffectLibrary::CountOfFreeGraphs() {
		return free_global_graphs.size() + free_property_graphs.size() + free_base_graphs.size() + free_variation_graphs.size() + free_overtime_graphs.size();
	}

	DataTypesDictionary::DataTypesDictionary() {
		eff.Insert("name", tfxString);
		eff.Insert("image_index", tfxUint);
		eff.Insert("image_handle_x", tfxFloat);
		eff.Insert("image_handle_y", tfxFloat);
		eff.Insert("spawn_amount", tfxUint);
		eff.Insert("blend_mode", tfxSInt);
		eff.Insert("image_start_frame", tfxFloat);
		eff.Insert("image_end_frame", tfxFloat);

		eff.Insert("emission_type", tfxSInt);
		eff.Insert("emission_direction", tfxSInt);
		eff.Insert("grid_rows", tfxFloat);
		eff.Insert("grid_columns", tfxFloat);
		eff.Insert("loop_length", tfxFloat);
		eff.Insert("emitter_handle_x", tfxFloat);
		eff.Insert("emitter_handle_y", tfxFloat);
		eff.Insert("end_behaviour", tfxSInt);
		eff.Insert("angle_setting", tfxSInt);
		eff.Insert("angle_offset", tfxFloat);
		eff.Insert("multiply_blend_factor", tfxFloat);

		eff.Insert("random_color", tfxBool);
		eff.Insert("relative_position", tfxBool);
		eff.Insert("relative_angle", tfxBool);
		eff.Insert("image_handle_auto_center", tfxBool);
		eff.Insert("single", tfxBool);
		eff.Insert("one_shot", tfxBool);
		eff.Insert("spawn_on_grid", tfxBool);
		eff.Insert("grid_spawn_clockwise", tfxBool);
		eff.Insert("fill_area", tfxBool);
		eff.Insert("emitter_handle_auto_center", tfxBool);
		eff.Insert("edge_traversal", tfxBool);
		eff.Insert("image_reverse_animation", tfxBool);
		eff.Insert("image_play_once", tfxBool);
		eff.Insert("image_animate", tfxBool);
		eff.Insert("image_random_start_frame", tfxBool);
		eff.Insert("global_uniform_size", tfxBool);
		eff.Insert("base_uniform_size", tfxBool);
		eff.Insert("lifetime_uniform_size", tfxBool);

		eff.Insert("frames", tfxUint);
		eff.Insert("current_frame", tfxUint);
		eff.Insert("frame_offset", tfxUint);
		eff.Insert("layer", tfxUint);
		eff.Insert("position_x", tfxFloat);
		eff.Insert("position_y", tfxFloat);
		eff.Insert("frame_width", tfxFloat);
		eff.Insert("frame_height", tfxFloat);
		eff.Insert("loop", tfxBool);
		eff.Insert("seamless", tfxBool);
		eff.Insert("seed", tfxUint);
		eff.Insert("zoom", tfxFloat);
		eff.Insert("scale", tfxFloat);
		eff.Insert("color_option", tfxSInt);
		eff.Insert("export_option", tfxSInt);
		eff.Insert("export_with_transparency", tfxBool);

		//Editor config, move this to the editor
		eff.Insert("only_play_selected_emitter", tfxBool);
		eff.Insert("load_examples", tfxBool);
		eff.Insert("load_last_file", tfxBool);
		eff.Insert("load_last_file_path", tfxString);
		eff.Insert("recent1", tfxString);
		eff.Insert("recent2", tfxString);
		eff.Insert("recent3", tfxString);
		eff.Insert("recent4", tfxString);
		eff.Insert("background_color_red", tfxFloat);
		eff.Insert("background_color_green", tfxFloat);
		eff.Insert("background_color_blue", tfxFloat);
		eff.Insert("use_checker_background", tfxBool);
		eff.Insert("preview_zoom", tfxFloat);
		eff.Insert("updates_per_second", tfxFloat);
		eff.Insert("background_image", tfxString);
		eff.Insert("use_background_image", tfxBool);
		eff.Insert("background_image_scale_x", tfxFloat);
		eff.Insert("background_image_scale_y", tfxFloat);
		eff.Insert("background_image_offset_x", tfxFloat);
		eff.Insert("background_image_offset_y", tfxFloat);
		eff.Insert("autoplay_effect", tfxSInt);
		eff.Insert("sync_refresh_rate", tfxBool);
		eff.Insert("window_maximised", tfxBool);
		eff.Insert("window_width", tfxSInt);
		eff.Insert("window_height", tfxSInt);
		eff.Insert("show_emitter_positions", tfxBool);
		eff.Insert("dpi_factor", tfxFloat);
		eff.Insert("graph_lookup_mode", tfxSInt);
		eff.Insert("show_tool_tips", tfxBool);
		eff.Insert("preview_trail_mode", tfxBool);
	}

	int ValidateEffectLibrary(const char *filename) {
		mz_zip_archive zip_archive;
		memset(&zip_archive, 0, sizeof(zip_archive));
		int error = 0;

		bool status = mz_zip_reader_init_file(&zip_archive, filename, 0);
		if (!status)
		{
			error = -1;
		}

		if ((int)mz_zip_reader_get_num_files(&zip_archive) == 0)
			error = -1;

		size_t uncomp_size;
		std::stringstream d;
		void *data = mz_zip_reader_extract_file_to_heap(&zip_archive, "data.txt", &uncomp_size, 0);

		if (!data)
		{
			error = -1;
		}

		mz_free(data);
		d.clear();
		mz_zip_reader_end(&zip_archive);

		return error;
	}

	void AssignGraphData(EffectEmitter &effect, tfxvec<tfxText> &values) {
		if (values.size() > 0) {
			if (values[0] == "global_amount") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].amount.AddNode(n); }
			if (values[0] == "global_effect_angle") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].effect_angle.AddNode(n); }
			if (values[0] == "global_frame_rate") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].frame_rate.AddNode(n); }
			if (values[0] == "global_height") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].height.AddNode(n); }
			if (values[0] == "global_width") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].width.AddNode(n); }
			if (values[0] == "global_life") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].life.AddNode(n); }
			if (values[0] == "global_opacity") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].opacity.AddNode(n); }
			if (values[0] == "global_spin") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].spin.AddNode(n); }
			if (values[0] == "global_splatter") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].splatter.AddNode(n); }
			if (values[0] == "global_stretch") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].stretch.AddNode(n); }
			if (values[0] == "global_overal_scale") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].overal_scale.AddNode(n); }
			if (values[0] == "global_weight") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].weight.AddNode(n); }
			if (values[0] == "global_velocity") { AttributeNode n; AssignNodeData(n, values); effect.library->global_graphs[effect.global].velocity.AddNode(n); }

			if (values[0] == "base_arc_offset") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].arc_offset.AddNode(n); }
			if (values[0] == "base_arc_size") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].arc_size.AddNode(n); }
			if (values[0] == "base_emission_angle") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emission_angle.AddNode(n); }
			if (values[0] == "base_emission_range") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emission_range.AddNode(n); }
			if (values[0] == "base_emitter_height") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emitter_height.AddNode(n); }
			if (values[0] == "base_emitter_width") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emitter_width.AddNode(n); }
			if (values[0] == "base_splatter") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].splatter.AddNode(n); }

			if (values[0] == "property_arc_offset") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].arc_offset.AddNode(n); }
			if (values[0] == "property_arc_size") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].arc_size.AddNode(n); }
			if (values[0] == "property_emitter_angle") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emitter_angle.AddNode(n); }
			if (values[0] == "property_emission_angle") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emission_angle.AddNode(n); }
			if (values[0] == "property_emission_range") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emission_range.AddNode(n); }
			if (values[0] == "property_emitter_height") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emitter_height.AddNode(n); }
			if (values[0] == "property_emitter_width") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].emitter_width.AddNode(n); }
			if (values[0] == "property_splatter") { AttributeNode n; AssignNodeData(n, values); effect.library->property_graphs[effect.property].splatter.AddNode(n); }

			if (values[0] == "base_amount") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].amount.AddNode(n); }
			if (values[0] == "base_life") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].life.AddNode(n); }
			if (values[0] == "base_height") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].height.AddNode(n); }
			if (values[0] == "base_width") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].width.AddNode(n); }
			if (values[0] == "base_spin") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].spin.AddNode(n); }
			if (values[0] == "base_velocity") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].velocity.AddNode(n); }
			if (values[0] == "base_weight") { AttributeNode n; AssignNodeData(n, values); effect.library->base_graphs[effect.base].weight.AddNode(n); }

			if (values[0] == "variation_amount") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].amount.AddNode(n); }
			if (values[0] == "variation_height") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].height.AddNode(n); }
			if (values[0] == "variation_width") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].width.AddNode(n); }
			if (values[0] == "variation_life") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].life.AddNode(n); }
			if (values[0] == "variation_velocity") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].velocity.AddNode(n); }
			if (values[0] == "variation_weight") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].weight.AddNode(n); }
			if (values[0] == "variation_spin") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].spin.AddNode(n); }
			if (values[0] == "variation_motion_randomness") { AttributeNode n; AssignNodeData(n, values); effect.library->variation_graphs[effect.variation].motion_randomness.AddNode(n); }

			if (values[0] == "overtime_red") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].red.AddNode(n); }
			if (values[0] == "overtime_green") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].green.AddNode(n); }
			if (values[0] == "overtime_blue") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].blue.AddNode(n); }
			if (values[0] == "overtime_opacity") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].opacity.AddNode(n); }
			if (values[0] == "overtime_intensity") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].intensity.AddNode(n); }
			if (values[0] == "overtime_frame_rate") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].frame_rate.AddNode(n); }
			if (values[0] == "overtime_spin") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].spin.AddNode(n); }
			if (values[0] == "overtime_stretch") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].stretch.AddNode(n); }
			if (values[0] == "overtime_velocity") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].velocity.AddNode(n); }
			if (values[0] == "overtime_weight") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].weight.AddNode(n); }
			if (values[0] == "overtime_width") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].width.AddNode(n); }
			if (values[0] == "overtime_height") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].height.AddNode(n); }
			if (values[0] == "overtime_motion_randomness") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].motion_randomness.AddNode(n); }
			if (values[0] == "overtime_direction") { AttributeNode n; AssignNodeData(n, values); effect.library->overtime_graphs[effect.overtime].direction.AddNode(n); }
		}
	}

	void AssignNodeData(AttributeNode &n, tfxvec<tfxText> &values) {
		n.frame = (float)atof(values[1].c_str());
		n.value = (float)atof(values[2].c_str());
		n.is_curve = (bool)atoi(values[3].c_str());
		n.left.x = (float)atof(values[4].c_str());
		n.left.y = (float)atof(values[5].c_str());
		n.right.x = (float)atof(values[6].c_str());
		n.right.y = (float)atof(values[7].c_str());
		if (n.is_curve)
			n.curves_initialised = true;
	}

	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, uint32_t value) {
		if (field == "image_index")
			effect.properties.shape_index = value;
		if (field == "spawn_amount")
			effect.properties.spawn_amount = value;
		if (field == "frames")
			effect.library->animation_settings[effect.animation_settings].frames = value;
		if (field == "current_frame")
			effect.library->animation_settings[effect.animation_settings].current_frame = value;
		if (field == "seed")
			effect.library->animation_settings[effect.animation_settings].seed = value;
		if (field == "layer")
			effect.properties.layer = value;
	}
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, int value) {
		if (field == "emission_type")
			effect.properties.emission_type = (EmissionType)value;
		if (field == "emission_direction")
			effect.properties.emission_direction = (EmissionDirection)value;
		if (field == "blend_mode")
			effect.properties.blend_mode = (BlendMode)value;
		if (field == "angle_setting")
			effect.properties.angle_setting = (AngleSetting)value;
		if (field == "color_option")
			effect.library->animation_settings[effect.animation_settings].color_option = (ExportColorOptions)value;
		if (field == "export_option")
			effect.library->animation_settings[effect.animation_settings].export_option = (ExportOptions)value;
		if (field == "end_behaviour")
			effect.properties.end_behaviour = (LineTraversalEndBehaviour)value;
	}
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, tfxText &value) {
		if (field == "name")
			effect.name = value;
	}
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, float value) {
		if (field == "position_x")
			effect.library->animation_settings[effect.animation_settings].position.x = value;
		if (field == "position_y")
			effect.library->animation_settings[effect.animation_settings].position.y = value;
		if (field == "frame_width")
			effect.library->animation_settings[effect.animation_settings].frame_size.x = value;
		if (field == "frame_height")
			effect.library->animation_settings[effect.animation_settings].frame_size.y = value;
		if (field == "zoom")
			effect.library->animation_settings[effect.animation_settings].zoom = value;
		if (field == "scale")
			effect.library->animation_settings[effect.animation_settings].scale = value;
		if (field == "image_handle_x")
			effect.properties.image_handle.x = value;
		if (field == "image_handle_y")
			effect.properties.image_handle.y = value;
		if (field == "grid_rows")
			effect.properties.grid_points.x = value;
		if (field == "grid_columns")
			effect.properties.grid_points.y = value;
		if (field == "loop_length")
			effect.properties.loop_length = value;
		if (field == "emitter_handle_x")
			effect.properties.emitter_handle.x = value;
		if (field == "emitter_handle_y")
			effect.properties.emitter_handle.y = value;
		if (field == "image_start_frame")
			effect.properties.start_frame = value;
		if (field == "image_end_frame")
			effect.properties.end_frame = value;
		if (field == "angle_offset")
			effect.properties.angle_offset = value;
	}
	void AssignEffectorProperty(EffectEmitter &effect, tfxText &field, bool value) {
		if (field == "loop")
			effect.library->animation_settings[effect.animation_settings].loop = value;
		if (field == "seamless")
			effect.library->animation_settings[effect.animation_settings].seamless = value;
		if (field == "export_with_transparency")
			effect.library->animation_settings[effect.animation_settings].export_with_transparency = value;
		if (field == "random_color")
			if (value) effect.properties.flags |= tfxEmitterPropertyFlags_random_color; else effect.properties.flags &= ~tfxEmitterPropertyFlags_random_color;
		if (field == "relative_position")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_relative_position; else effect.properties.flags &= ~tfxEmitterPropertyFlags_relative_position;
		if (field == "relative_angle")
			if (value) effect.properties.flags |= tfxEmitterPropertyFlags_relative_angle; else effect.properties.flags &= ~tfxEmitterPropertyFlags_relative_angle;
		if (field == "image_handle_auto_center")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_image_handle_auto_center; else effect.properties.flags &= ~tfxEmitterPropertyFlags_image_handle_auto_center;
		if (field == "single")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_single; else effect.properties.flags &= ~tfxEmitterPropertyFlags_single;
		if (field == "one_shot")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_one_shot; else effect.properties.flags &= ~tfxEmitterPropertyFlags_one_shot;
		if (field == "spawn_on_grid")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_spawn_on_grid; else effect.properties.flags &= ~tfxEmitterPropertyFlags_spawn_on_grid;
		if (field == "grid_spawn_clockwise")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_grid_spawn_clockwise; else effect.properties.flags &= ~tfxEmitterPropertyFlags_grid_spawn_clockwise;
		if (field == "fill_area")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_fill_area; else effect.properties.flags &= ~tfxEmitterPropertyFlags_fill_area;
		if (field == "emitter_handle_auto_center")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_emitter_handle_auto_center; else effect.properties.flags &= ~tfxEmitterPropertyFlags_emitter_handle_auto_center;
		if (field == "edge_traversal")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_edge_traversal; else effect.properties.flags &= ~tfxEmitterPropertyFlags_edge_traversal;
		if (field == "image_reverse_animation")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_reverse_animation; else effect.properties.flags &= ~tfxEmitterPropertyFlags_reverse_animation;
		if (field == "image_play_once")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_play_once; else effect.properties.flags &= ~tfxEmitterPropertyFlags_play_once;
		if (field == "image_animate")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_animate; else effect.properties.flags &= ~tfxEmitterPropertyFlags_animate;
		if (field == "image_random_start_frame")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_random_start_frame; else effect.properties.flags &= ~tfxEmitterPropertyFlags_random_start_frame;
		if (field == "global_uniform_size")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_global_uniform_size; else effect.properties.flags &= ~tfxEmitterPropertyFlags_global_uniform_size;
		if (field == "base_uniform_size")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_base_uniform_size; else effect.properties.flags &= ~tfxEmitterPropertyFlags_base_uniform_size;
		if (field == "lifetime_uniform_size")
			if(value) effect.properties.flags |= tfxEmitterPropertyFlags_lifetime_uniform_size; else effect.properties.flags &= ~tfxEmitterPropertyFlags_lifetime_uniform_size;
	}

	void StreamProperties(EmitterProperties &property, std::stringstream &file) {

		file << "image_index" Del property.shape_index EndLine;
		file << "image_handle_x" Del property.image_handle.x EndLine;
		file << "image_handle_y" Del property.image_handle.y EndLine;
		file << "image_start_frame" Del property.start_frame EndLine;
		file << "image_end_frame" Del property.end_frame EndLine;
		file << "image_play_once" Del (property.flags & tfxEmitterPropertyFlags_play_once) EndLine;
		file << "image_reverse_animation" Del (property.flags & tfxEmitterPropertyFlags_reverse_animation) EndLine;
		file << "image_animate" Del (property.flags & tfxEmitterPropertyFlags_animate) EndLine;
		file << "image_random_start_frame" Del (property.flags & tfxEmitterPropertyFlags_random_start_frame) EndLine;
		file << "spawn_amount" Del property.spawn_amount EndLine;
		file << "blend_mode" Del property.blend_mode EndLine;
		file << "emission_type" Del property.emission_type EndLine;
		file << "emission_direction" Del property.emission_direction EndLine;
		file << "grid_rows" Del property.grid_points.x EndLine;
		file << "grid_columns" Del property.grid_points.y EndLine;
		file << "loop_length" Del property.loop_length EndLine;
		file << "emitter_handle_x" Del property.emitter_handle.x EndLine;
		file << "emitter_handle_y" Del property.emitter_handle.y EndLine;
		file << "end_behaviour" Del property.end_behaviour EndLine;
		file << "random_color" Del (property.flags & tfxEmitterPropertyFlags_random_color) EndLine;
		file << "relative_position" Del (property.flags & tfxEmitterPropertyFlags_relative_position) EndLine;
		file << "relative_angle" Del (property.flags & tfxEmitterPropertyFlags_relative_angle) EndLine;
		file << "image_handle_auto_center" Del (property.flags & tfxEmitterPropertyFlags_image_handle_auto_center) EndLine;
		file << "single" Del (property.flags & tfxEmitterPropertyFlags_single) EndLine;
		file << "one_shot" Del (property.flags & tfxEmitterPropertyFlags_one_shot) EndLine;
		file << "spawn_on_grid" Del (property.flags & tfxEmitterPropertyFlags_spawn_on_grid) EndLine;
		file << "grid_spawn_clockwise" Del (property.flags & tfxEmitterPropertyFlags_grid_spawn_clockwise) EndLine;
		file << "fill_area" Del (property.flags & tfxEmitterPropertyFlags_fill_area) EndLine;
		file << "emitter_handle_auto_center" Del (property.flags & tfxEmitterPropertyFlags_emitter_handle_auto_center) EndLine;
		file << "edge_traversal" Del (property.flags & tfxEmitterPropertyFlags_edge_traversal) EndLine;
		file << "angle_setting" Del property.angle_setting EndLine;
		file << "angle_offset" Del property.angle_offset EndLine;
		file << "global_uniform_size" Del (property.flags & tfxEmitterPropertyFlags_global_uniform_size) EndLine;
		file << "base_uniform_size" Del (property.flags & tfxEmitterPropertyFlags_base_uniform_size) EndLine;
		file << "lifetime_uniform_size" Del (property.flags & tfxEmitterPropertyFlags_lifetime_uniform_size) EndLine;
		file << "layer" Del property.layer EndLine;

	}

	void StreamGraph(const char * name, Graph &graph, std::stringstream &file) {

		for (auto &n : graph.nodes) {
			file << name Com n.frame Com n.value Com n.is_curve Com n.left.x Com n.left.y Com n.right.x Com n.right.y EndLine;
		}

	}

	Random::Random() {
		ReSeed();
	}

	void Random::ReSeed() {
		seeds[0] = (uint64_t)Millisecs();
		seeds[1] = (uint64_t)Millisecs() * 2;
	}

	void Random::ReSeed(uint64_t seed1, uint64_t seed2) {
		seeds[0] = seed1;
		seeds[1] = seed2;
	}

	double Random::Millisecs() {
		auto now = Clock::now().time_since_epoch();
		auto m = std::chrono::duration_cast<std::chrono::milliseconds>(now).count();
		return double(m);
	}

	static bool CompareNodes(AttributeNode &left, AttributeNode &right) {
		return left.frame < right.frame;
	}

	bool Graph::IsOvertimeGraph() {
		return type >= tfxOvertime_velocity && type != tfxOvertime_velocity_adjuster;
	}

	bool Graph::IsGlobalGraph() {
		return type >= tfxGlobal_life && type <= tfxGlobal_splatter;
	}

	bool Graph::IsAngleGraph() {
		return (type == tfxGlobal_effect_angle || type == tfxProperty_emission_angle || type == tfxProperty_emission_range || type == tfxProperty_emitter_angle ||
			type == tfxProperty_arc_offset || type == tfxProperty_arc_size || type == tfxBase_spin || type == tfxVariation_spin || type == tfxOvertime_direction);
	}

	float AttributeNode::GetX() {
		return frame;
	}
	float AttributeNode::GetY() {
		return value;
	}

	bool SetNode(Graph &graph, AttributeNode &node, float _frame, float _value, bool _is_curve, float _c0x, float _c0y, float _c1x, float _c1y) {
		node.frame = _frame;
		node.value = _value;
		node.is_curve = _is_curve;
		node.left.x = _c0x;
		node.left.y = _c0y;
		node.right.x = _c1x;
		node.right.y = _c1y;
		if (graph.nodes[0] == node) {
			node.frame = graph.min.x;
			ClampNode(graph, node);
		}
		else {
			ClampNode(graph, node);
		}

		if (graph.Sort()) {
			graph.ReIndex();
			return true;
		}

		return false;
	}

	bool SetNodeFrame(Graph &graph, AttributeNode &node, float &_frame) {
		node.frame = _frame;
		if (graph.nodes[0] == node) {
			node.frame = graph.min.x;
			ClampNode(graph, node);
		}
		else {
			ClampNode(graph, node);
		}
		_frame = node.frame;

		if (graph.Sort()) {
			graph.ReIndex();
			return true;
		}

		return false;

	}

	bool SetNodeValue(Graph &graph, AttributeNode &node, float &_value) {
		node.value = _value;
		ClampNode(graph, node);
		_value = node.value;

		return false;
	}

	bool SetNode(Graph &graph, AttributeNode &node, float &frame, float &value) {
		float old_frame = node.frame;
		float old_value = node.value;

		node.frame = frame;
		node.value = value;

		if (graph.nodes[0] == node) {
			node.frame = graph.min.x;
			ClampNode(graph, node);
		}
		else {
			ClampNode(graph, node);
		}

		if (node.curves_initialised) {
			node.left.y += node.value - old_value;
			node.left.x += node.frame - old_frame;
			node.right.y += node.value - old_value;
			node.right.x += node.frame - old_frame;
		}

		frame = node.frame;
		value = node.value;

		if (graph.Sort()) {
			graph.ReIndex();
			return true;
		}

		return false;
	}

	void SetCurve(Graph &graph, AttributeNode &node, bool is_left_curve, float &frame, float &value) {
		if (is_left_curve) {
			node.left.x = frame;
			node.left.y = value;
			if (node.left.x > node.frame)
				node.left.x = node.frame;
			else
				ClampCurve(graph, node.left, node);
			frame = node.left.x;
			value = node.left.y;
		}
		else {
			node.right.x = frame;
			node.right.y = value;
			if (node.right.x < node.frame)
				node.right.x = node.frame;
			else
				ClampCurve(graph, node.right, node);
			frame = node.right.x;
			value = node.right.y;
		}
	}

	bool MoveNode(Graph &graph, AttributeNode &node, float frame, float value, bool sort) {
		float old_frame = node.frame;
		float old_value = node.value;

		node.frame += frame;
		node.value += value;

		if (graph.nodes[0] == node) {
			node.frame = graph.min.x;
			ClampNode(graph, node);
		}
		else {
			ClampNode(graph, node);
		}

		if (node.curves_initialised) {
			node.left.y += node.value - old_value;
			node.left.x += node.frame - old_frame;
			node.right.y += node.value - old_value;
			node.right.x += node.frame - old_frame;
		}

		if (sort) {
			if (graph.Sort()) {
				graph.ReIndex();
				return true;
			}
		}

		return false;
	}

	void ClampNode(Graph &graph, AttributeNode &node) {
		if (node.value < graph.min.y) node.value = graph.min.y;
		if (node.frame < graph.min.x) node.frame = graph.min.x;
		if (node.value > graph.max.y) node.value = graph.max.y;
		if (node.frame > graph.max.x) node.frame = graph.max.x;
	}

	void ClampCurve(Graph &graph, Point &p, AttributeNode &node) {
		if (p.y < graph.min.y) p.y = graph.min.y;
		if (p.x < graph.min.x) p.x = graph.min.x;
		//if (p.y > graph.max.y) p.y = graph.max.y;
		if (p.x > graph.max.x) p.x = graph.max.x;

		AttributeNode *next = graph.GetNextNode(node);
		if (next) {
			if (p.x > next->frame) p.x = next->frame;
		}

		AttributeNode *prev = graph.GetPrevNode(node);
		if (prev) {
			if (p.x < prev->frame) p.x = prev->frame;
		}
	}

	Graph::Graph() {
		min.x = 0.f;
		min.y = 0.f;
		max.x = 1000.f;
		max.y = 1000.f;

		effector = nullptr;
	}

	Graph::~Graph() {
		Free();
	}

	const float BEZIER_ACCURACY = 0.01f;

	float GetDistance(float fromx, float fromy, float tox, float toy) {

		float w = tox - fromx;
		float h = toy - fromy;

		return std::sqrt(w * w + h * h);

	}

	float GetVectorAngle(float x, float y) {
		return std::atan2f(x, -y);
	}

	Point GetQuadBezier(Point p0, Point p1, Point p2, float t, float ymin, float ymax, bool clamp) {
		Point b;
		b.x = std::powf(1.f - t, 2.f) * p0.x + 2.f * t * (1.f - t) * p1.x + std::powf(t, 2.f) * p2.x;
		b.y = std::powf(1.f - t, 2.f) * p0.y + 2.f * t * (1.f - t) * p1.y + std::powf(t, 2.f) * p2.y;
		if (b.x < p0.x) b.x = p0.x;
		if (b.x > p2.x) b.x = p2.x;
		if (clamp) {
			if (b.y < ymin) b.y = ymin;
			if (b.y > ymax) b.y = ymax;
		}
		return b;
	}

	Point GetCubicBezier(Point p0, Point p1, Point p2, Point p3, float t, float ymin, float ymax, bool clamp) {
		Point b;
		b.x = std::powf(1.f - t, 3.f) * p0.x + 3.f * t * std::powf(1.f - t, 2.f) * p1.x + 3.f * std::powf(t, 2.f) * (1.f - t) * p2.x + std::powf(t, 3.f) * p3.x;
		b.y = std::powf(1.f - t, 3.f) * p0.y + 3.f * t * std::powf(1.f - t, 2.f) * p1.y + 3.f * std::powf(t, 2.f) * (1.f - t) * p2.y + std::powf(t, 3.f) * p3.y;
		if (b.x < p0.x) b.x = p0.x;
		if (b.x > p3.x) b.x = p3.x;
		if (clamp) {
			if (b.y < ymin) b.y = ymin;
			if (b.y > ymax) b.y = ymax;
		}
		return b;
	}

	float GetBezierValue(const AttributeNode *lastec, const AttributeNode &a, float t, float ymin, float ymax) {
		if (lastec) {
			if (a.is_curve) {
				if (lastec->is_curve) {
					Point p0(lastec->frame, lastec->value);
					Point p1(lastec->right.x, lastec->right.y);
					Point p2(a.left.x, a.left.y);
					Point p3(a.frame, a.value);
					Point value = GetCubicBezier(p0, p1, p2, p3, t, ymin, ymax);
					return value.y;
				}
				else {
					Point p0(lastec->frame, lastec->value);
					Point p1(a.left.x, a.left.y);
					Point p2(a.frame, a.value);
					Point value = GetQuadBezier(p0, p1, p2, t, ymin, ymax);
					return value.y;
				}
			}
			else if (lastec->is_curve) {
				Point p0(lastec->frame, lastec->value);
				Point p1(lastec->right.x, lastec->right.y);
				Point p2(a.frame, a.value);
				Point value = GetQuadBezier(p0, p1, p2, t, ymin, ymax);
				return value.y;
			}
		}
		else {
			return 0;
		}

		return 0;
	}

	AttributeNode* Graph::AddNode(float _frame, float _value, bool _is_curve, float _c0x, float _c0y, float _c1x, float _c1y) {
		AttributeNode node;

		if (nodes.size())
			node.frame = _frame;
		else
			node.frame = 0.f;

		node.value = _value;
		node.is_curve = _is_curve;
		node.left.x = _c0x;
		node.left.y = _c0y;
		node.right.x = _c1x;
		node.right.y = _c1y;
		ClampNode(*this, node);
		nodes.push_back(node);
		Sort();

		ReIndex();
		return &nodes.back();
	}

	void Graph::AddNode(AttributeNode &node) {
		for (auto &n : nodes) {
			if (n.frame == node.frame)
				return;
		}
		nodes.push_back(node);
		Sort();
		ReIndex();
	}

	AttributeNode* Graph::AddCoordNode(float _frame, float _value) {
		AttributeNode node;

		if (nodes.size())
			node.frame = _frame;
		else
			node.frame = 0.f;

		node.value = _value;
		node.is_curve = false;
		node.left.x = 0.f;
		node.left.y = 0.f;
		node.right.x = 0.f;
		node.right.y = 0.f;
		ClampNode(*this, node);
		AttributeNode &n = nodes.push_back(node);
		if (Sort()) {
			ReIndex();
			return nodes.find(n);
		}

		ReIndex();
		return &n;
	}

	AttributeNode* Graph::InsertCoordNode(float _frame, float _value) {
		AttributeNode node;

		if (nodes.size())
			node.frame = _frame;
		else
			node.frame = 0.f;

		node.value = _value;
		node.is_curve = false;
		node.left.x = 0.f;
		node.left.y = 0.f;
		node.right.x = 0.f;
		node.right.y = 0.f;
		ClampNode(*this, node);

		if (nodes.size() > 1) {
			AttributeNode *last_node = nullptr;
			for (auto *n = nodes.begin() + 1; n != nodes.end(); ++n) {
				if (node.frame < n->frame)
					last_node = n;
				else
					break;
			}

			if (last_node) {
				AttributeNode *r_value = nodes.insert(last_node, node);
				ReIndex();
				return r_value;
			}
		}

		AttributeNode *r_value = &nodes.push_back(node);
		ReIndex();
		return r_value;
	}

	AttributeNode* Graph::InsertNode(float _frame, float _value) {
		AttributeNode node;

		if (nodes.size())
			node.frame = _frame;
		else
			node.frame = 0.f;

		node.value = _value;
		node.is_curve = false;
		node.left.x = 0.f;
		node.left.y = 0.f;
		node.right.x = 0.f;
		node.right.y = 0.f;
		ClampNode(*this, node);

		if (nodes.size() > 1) {
			AttributeNode *last_node = nullptr;
			for (auto *n = nodes.begin() + 1; n != nodes.end(); ++n) {
				if (node.frame < n->frame)
					last_node = n;
				else
					break;
			}

			if (last_node) {
				AttributeNode *r_value = nodes.insert(last_node, node);
				ReIndex();
				return r_value;
			}
		}

		AttributeNode *r_value = &nodes.push_back(node);
		ReIndex();
		return r_value;
	}

	void Graph::SetNode(uint32_t i, float _frame, float _value, bool _is_curve, float _c0x, float _c0y, float _c1x, float _c1y) {
		if (!nodes.empty() && i < nodes.size()) {
			nodes[i].frame = _frame;
			nodes[i].value = _value;
			nodes[i].is_curve = _is_curve;
			nodes[i].left.x = _c0x;
			nodes[i].left.y = _c0y;
			nodes[i].right.x = _c1x;
			nodes[i].right.y = _c1y;
			if (Sort())
				ReIndex();
		}
	}

	tfxvec<AttributeNode>& Graph::Nodes() {
		return nodes;
	}

	float Graph::GetValue(float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		for (auto &a : nodes) {
			if (age < a.frame) {
				p = (age - lastf) / (a.frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, min.y, max.y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = a.frame - 1;
			lastec = &a;
		}
		return lastv;

	}

	AttributeNode *Graph::GetNextNode(AttributeNode &node) {
		if (node.index < nodes.size() - 1) {
			return &nodes[node.index + 1];
		}

		return nullptr;
	}

	AttributeNode *Graph::GetPrevNode(AttributeNode &node) {
		if (node.index > 0) {
			return &nodes[node.index - 1];
		}

		return nullptr;
	}

	AttributeNode *Graph::GetLastNode() {
		return &nodes.back();
	}

	float Graph::GetRandomValue(float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		for (auto &a : nodes) {
			if (age < a.frame) {
				p = (age - lastf) / (a.frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, min.y, max.y);
				if (bezier_value) {
					return random_generation.Range(bezier_value);
				}
				else {
					return random_generation.Range(lastv - p * (lastv - a.value));
				}
			}
			lastv = a.value;
			lastf = a.frame - 1;
			lastec = &a;
		}
		return random_generation.Range(lastv);

	}

	float Graph::GetValue(float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		for (auto &a : nodes) {
			float frame = a.frame * life;
			if (age < frame) {
				p = (age - lastf) / (frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, min.y, max.y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = frame - 1;
			lastec = &a;
		}
		return lastv;
	}

	float Graph::GetFirstValue() {
		if (nodes.size())
			return nodes.front().value;
		return 0.f;
	}

	float* Graph::LinkFirstValue() {
		if (nodes.size())
			return &nodes.front().value;
		return nullptr;
	}

	float Graph::GetLastValue() {
		if (nodes.size())
			return nodes.back().value;

		return 0.f;
	}

	float Graph::GetMaxValue() {
		if (nodes.size()) {
			float value = tfxMIN_FLOAT;
			for (auto &n : nodes) {
				if (value < n.value)
					value = n.value;
			}
			return value;
		}
		return 0.f;
	}

	float Graph::GetMinValue() {
		if (nodes.size()) {
			float value = tfxMAX_FLOAT;
			for (auto &n : nodes) {
				if (value > n.value)
					value = n.value;
			}
			return value;
		}
		return 0.f;
	}

	float Graph::GetLastFrame() {
		if (nodes.size())
			return nodes.back().frame;

		return 0.f;
	}

	AttributeNode* Graph::FindNode(const AttributeNode &n) {
		return nodes.find(n);
	}

	void Graph::ValidateCurves() {
		unsigned int index = 0;
		unsigned int last_index = nodes.size() - 1;
		for(auto &n : nodes) {
			if (n.is_curve) {
				if (index < last_index) {
					if (nodes[index + 1].frame < n.right.x)
						n.right.x = nodes[index + 1].frame;
				}
				if (index > 0) {
					if (nodes[index - 1].frame > n.left.x)
						n.left.x = nodes[index - 1].frame;
				}
				if (n.left.x > n.frame)
					n.left.x = n.frame;
				if (n.right.x < n.frame)
					n.right.x = n.frame;
			}
			index++;
		}
	}

	void Graph::DeleteNode(const AttributeNode &n) {
		nodes.erase(&n);
	}

	void Graph::Reset(float v, GraphPreset preset, bool add_node) {
		nodes.clear();
		if (add_node)
			AddNode(0.f, v);
		switch (preset) {
		case GraphPreset::tfxGlobalPercentPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 20.f };
			break;
		case GraphPreset::tfxGlobalOpacityPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 1.f };
			break;
		case GraphPreset::tfxGlobalPercentPresetSigned:
			min = { 0.f, -20.f }; max = { tfxMAX_FRAME, 20.f };
			break;
		case GraphPreset::tfxAnglePreset:
			min = { 0.f, -1080.f }; max = { tfxMAX_FRAME, 1080.f };
			break;
		case GraphPreset::tfxArcPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 360.f };
			break;
		case GraphPreset::tfxEmissionRangePreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 360.f };
			break;
		case GraphPreset::tfxDimensionsPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 4000.f };
			break;
		case GraphPreset::tfxLifePreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 100000.f };
			break;
		case GraphPreset::tfxAmountPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 5000.f };
			break;
		case GraphPreset::tfxVelocityPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 10000.f };
			break;
		case GraphPreset::tfxWeightPreset:
			min = { 0.f, -2500.f }; max = { tfxMAX_FRAME, 2500.f };
			break;
		case GraphPreset::tfxWeightVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 2500.f };
			break;
		case GraphPreset::tfxSpinPreset:
			min = { 0.f, -2000.f }; max = { tfxMAX_FRAME, 2000.f };
			break;
		case GraphPreset::tfxSpinVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 2000.f };
			break;
		case GraphPreset::tfxDirectionVariationPreset:
			min = { 0.f, 0.f }; max = { tfxMAX_FRAME, 22.5f };
			break;
		case GraphPreset::tfxWeightOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 20.f };
			break;
		case GraphPreset::tfxDirectionOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 4320.f };
			break;
		case GraphPreset::tfxSpinOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 20.f };
			break;
		case GraphPreset::tfxVelocityOvertimePreset:
			min = { 0.f, -20.f }; max = { 1.f, 20.f };
			break;
		case GraphPreset::tfxPercentOvertime:
			min = { 0.f, 0.f }; max = { 1.f, 20.f };
			break;
		case GraphPreset::tfxFrameratePreset:
			min = { 0.f, 0.f }; max = { 1.f, 200.f };
			break;
		case GraphPreset::tfxOpacityOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 1.f };
			break;
		case GraphPreset::tfxColorPreset:
			min = { 0.f, 0.f }; max = { 1.f, 255.f };
			break;
		case GraphPreset::tfxIntensityOvertimePreset:
			min = { 0.f, 0.f }; max = { 1.f, 5.f };
			break;
		}

		graph_preset = preset;
	}

	tfxVec4 GetMinMaxGraphValues(GraphPreset preset) {
		tfxVec4 mm;
		switch (preset) {
		case GraphPreset::tfxGlobalPercentPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 20.f };
			break;
		case GraphPreset::tfxGlobalOpacityPreset:
			mm = { 0.f, 0.f , tfxMAX_FRAME, 1.f };
			break;
		case GraphPreset::tfxGlobalPercentPresetSigned:
			mm = { 0.f, -20.f, tfxMAX_FRAME, 20.f };
			break;
		case GraphPreset::tfxAnglePreset:
			mm = { 0.f, -1080.f, tfxMAX_FRAME, 1080.f };
			break;
		case GraphPreset::tfxArcPreset:
			mm = { 0.f, 0.f , tfxMAX_FRAME, 360.f };
			break;
		case GraphPreset::tfxEmissionRangePreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 360.f };
			break;
		case GraphPreset::tfxDimensionsPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 4000.f };
			break;
		case GraphPreset::tfxLifePreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 100000.f };
			break;
		case GraphPreset::tfxAmountPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 5000.f };
			break;
		case GraphPreset::tfxVelocityPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 10000.f };
			break;
		case GraphPreset::tfxWeightPreset:
			mm = { 0.f, -2500.f, tfxMAX_FRAME, 2500.f };
			break;
		case GraphPreset::tfxWeightVariationPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 2500.f };
			break;
		case GraphPreset::tfxSpinPreset:
			mm = { 0.f, -2000.f, tfxMAX_FRAME, 2000.f };
			break;
		case GraphPreset::tfxSpinVariationPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 2000.f };
			break;
		case GraphPreset::tfxDirectionVariationPreset:
			mm = { 0.f, 0.f, tfxMAX_FRAME, 22.5f };
			break;
		case GraphPreset::tfxWeightOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 20.f };
			break;
		case GraphPreset::tfxDirectionOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 4320.f };
			break;
		case GraphPreset::tfxSpinOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 20.f };
			break;
		case GraphPreset::tfxVelocityOvertimePreset:
			mm = { 0.f, -20.f, 1.f, 20.f };
			break;
		case GraphPreset::tfxPercentOvertime:
			mm = { 0.f, 0.f, 1.f, 20.f };
			break;
		case GraphPreset::tfxFrameratePreset:
			mm = { 0.f, 0.f, 1.f, 200.f };
			break;
		case GraphPreset::tfxOpacityOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 1.f };
			break;
		case GraphPreset::tfxColorPreset:
			mm = { 0.f, 0.f, 1.f, 255.f };
			break;
		case GraphPreset::tfxIntensityOvertimePreset:
			mm = { 0.f, 0.f, 1.f, 5.f };
			break;
		}
		
		return mm;
	}

	void Graph::DragValues(GraphPreset preset, float &frame, float &value) {
		switch (preset) {
		case GraphPreset::tfxOpacityOvertimePreset:
		case GraphPreset::tfxGlobalPercentPreset:
		case GraphPreset::tfxIntensityOvertimePreset:
			frame = 0.001f;
			value = 0.001f;
			break;
		case GraphPreset::tfxDirectionOvertimePreset:
			frame = 0.001f;
			value = 0.1f;
			break;
		case GraphPreset::tfxLifePreset:
			frame = 5;
			value = 5;
			break;
		case GraphPreset::tfxAnglePreset:
		case GraphPreset::tfxArcPreset:
		case GraphPreset::tfxEmissionRangePreset:
			frame = 5;
			value = 0.1f;
			break;
		case GraphPreset::tfxDimensionsPreset:
		case GraphPreset::tfxAmountPreset:
		case GraphPreset::tfxVelocityPreset:
		case GraphPreset::tfxWeightPreset:
		case GraphPreset::tfxWeightVariationPreset:
		case GraphPreset::tfxSpinPreset:
		case GraphPreset::tfxSpinVariationPreset:
		case GraphPreset::tfxFrameratePreset:
			frame = 5.f;
			value = 1.f;
			break;
		case GraphPreset::tfxWeightOvertimePreset:
		case GraphPreset::tfxVelocityOvertimePreset:
		case GraphPreset::tfxSpinOvertimePreset:
		case GraphPreset::tfxDirectionVariationPreset:
			frame = 0.001f;
			value = 0.01f;
			break;
		case GraphPreset::tfxColorPreset:
			frame = 0.001f;
			value = 1.f;
			break;
		default:
			frame = 1;
			value = 0.1f;
			break;
		}
	}

	void Graph::Clear() {
		nodes.clear();
	}

	void Graph::Free() {
		nodes.free_all();
	}

	void Graph::Copy(Graph &to) {
		to.nodes.reserve(nodes.size());
		std::copy(nodes.begin(), nodes.end(), to.nodes.begin());
		to.nodes.current_size = nodes.current_size;
	}

	bool Graph::Sort() {
		if (!std::is_sorted(nodes.begin(), nodes.end(), CompareNodes)) {
			std::sort(nodes.begin(), nodes.end(), CompareNodes);
			return true;
		}
		return false;
	}

	void Graph::ReIndex() {
		unsigned int i = 0;
		for (auto &a : nodes) {
			a.index = i++;
		}
	}

	tfxVec2 Graph::GetInitialZoom() {
		switch (graph_preset) {
		case GraphPreset::tfxOpacityOvertimePreset:
			return tfxVec2(.0017f, 0.00275f);
		case GraphPreset::tfxGlobalPercentPreset:
			return tfxVec2(3.5f, 0.005f);
			break;
		case GraphPreset::tfxGlobalPercentPresetSigned:
			return tfxVec2(3.5f, 0.006f);
			break;
		case GraphPreset::tfxGlobalOpacityPreset:
			return tfxVec2(3.5f, 0.003f);
			break;
		case GraphPreset::tfxLifePreset:
			return tfxVec2(3.5f, 3.5f);
			break;
		case GraphPreset::tfxAnglePreset:
			return tfxVec2(3.5f, 1.f);
			break;
		case GraphPreset::tfxArcPreset:
			return tfxVec2(3.5f, 1.f);
			break;
		case GraphPreset::tfxEmissionRangePreset:
			return tfxVec2(3.5f, .5f);
			break;
		case GraphPreset::tfxAmountPreset:
			return tfxVec2(3.5f, 1.25f);
			break;
		case GraphPreset::tfxFrameratePreset:
			return tfxVec2(0.0017f, .5f);
			break;
		case GraphPreset::tfxDimensionsPreset:
		case GraphPreset::tfxVelocityPreset:
		case GraphPreset::tfxWeightPreset:
		case GraphPreset::tfxWeightVariationPreset:
		case GraphPreset::tfxSpinPreset:
		case GraphPreset::tfxSpinVariationPreset:
			return tfxVec2(3.5f, 2.5f);
			break;
		case GraphPreset::tfxDirectionOvertimePreset:
			return tfxVec2(0.0017f, 1.f);
			break;
		case GraphPreset::tfxWeightOvertimePreset:
		case GraphPreset::tfxVelocityOvertimePreset:
		case GraphPreset::tfxSpinOvertimePreset:
		case GraphPreset::tfxDirectionVariationPreset:
		case GraphPreset::tfxPercentOvertime:
			return tfxVec2(0.0017f, 0.0035f);
			break;
		case GraphPreset::tfxIntensityOvertimePreset:
			return tfxVec2(0.0017f, 0.01115f);
			break;
		case GraphPreset::tfxColorPreset:
			break;
		default:
			return tfxVec2(0.1f, 0.1f);
			break;
		}

		return tfxVec2(0.1f, 0.1f);
	}

	void CompileGraph(Graph &graph) {
		float last_frame = graph.GetLastFrame();
		graph.lookup.last_frame = unsigned int(last_frame / tfxLOOKUP_FREQUENCY);
		graph.lookup.values.clear();
		if (graph.lookup.last_frame) {
			graph.lookup.values.resize(graph.lookup.last_frame + 1);
			for (unsigned int f = 0; f != graph.lookup.last_frame + 1; ++f) {
				graph.lookup.values[f] = graph.GetValue((float)f * tfxLOOKUP_FREQUENCY);
			}
			graph.lookup.values[graph.lookup.last_frame] = graph.GetLastValue();
		}
		else {
			graph.lookup.values.push_back(graph.GetFirstValue());
		}
	}

	void CompileGraphOvertime(Graph &graph) {
		graph.lookup.last_frame = unsigned int(graph.lookup.life / tfxLOOKUP_FREQUENCY_OVERTIME);
		graph.lookup.values.clear();
		if (graph.lookup.last_frame) {
			graph.lookup.values.resize(graph.lookup.last_frame + 1);
			for (unsigned int f = 0; f != graph.lookup.last_frame + 1; ++f) {
				graph.lookup.values[f] = graph.GetValue((float)f * tfxLOOKUP_FREQUENCY_OVERTIME, graph.lookup.life);
			}
			graph.lookup.values[graph.lookup.last_frame] = graph.GetLastValue();
		}
		else {
			graph.lookup.values.push_back(graph.GetFirstValue());
		}
	}

	float LookupFastOvertime(Graph &graph, float age, float lifetime) {
		float frame = 0;
		if (lifetime)
			frame = (age / lifetime * graph.lookup.life) / tfxLOOKUP_FREQUENCY_OVERTIME;
		if (frame < graph.lookup.last_frame)
			return graph.lookup.values[(unsigned int)frame];
		return graph.lookup.values[graph.lookup.last_frame];
	}

	float LookupFast(Graph &graph, float frame) {
		if ((unsigned int)frame < graph.lookup.last_frame)
			return graph.lookup.values[(unsigned int)frame];
		return graph.lookup.values[graph.lookup.last_frame];
	}

	float LookupPreciseOvertime(Graph &graph, float age, float life) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		for (auto &a : graph.nodes) {
			float frame = a.frame * life;
			if (age < frame) {
				p = (age - lastf) / (frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, graph.min.y, graph.max.y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = frame - 1;
			lastec = &a;
		}
		return lastv;
	}

	float LookupPrecise(Graph &graph, float age) {
		float lastv = 0;
		float lastf = 0;
		float p = 0;
		AttributeNode *lastec = nullptr;
		for (auto &a : graph.nodes) {
			if (age < a.frame) {
				p = (age - lastf) / (a.frame - lastf);
				float bezier_value = GetBezierValue(lastec, a, p, graph.min.y, graph.max.y);
				if (bezier_value) {
					return bezier_value;
				}
				else {
					return lastv - p * (lastv - a.value);
				}
			}
			lastv = a.value;
			lastf = a.frame - 1;
			lastec = &a;
		}
		return lastv;
	}

	float GetRandomFast(Graph &graph, float frame) {
		float value = 0;
		if ((unsigned int)frame < graph.lookup.last_frame)
			value = graph.lookup.values[(unsigned int)frame];
		value = graph.lookup.values[graph.lookup.last_frame];
		return random_generation.Range(value);
	}

	float GetRandomPrecise(Graph &graph, float frame) {
		return graph.GetRandomValue(frame);
	}

	float GetMaxLife(EffectEmitter &e) {
		Graph &life = *e.GetGraphByType(tfxBase_life);
		Graph &life_variation = *e.GetGraphByType(tfxVariation_life);
		float templife = 0;
		float max_life = 0;
		float life_last_frame = life.GetLastFrame();
		float life_variation_last_frame = life_variation.GetLastFrame();
		float global_adjust = 1.f;
		if (life_last_frame + life_variation_last_frame > 0) {
			for (float f = 0; f < std::fmaxf(life_last_frame, life_variation_last_frame); ++f) {
				if (e.parent)
					global_adjust = e.parent->GetGraphByType(tfxGlobal_life)->GetValue(f);
				templife = life.GetValue(f) + life_variation.GetValue(f);
				templife *= global_adjust;
				if (max_life < templife)
					max_life = templife;
			}
		}
		else {
			max_life = life.GetFirstValue() + life_variation.GetFirstValue();
		}

		return max_life;
	}

	bool IsOvertimeGraph(GraphType type) {
		return type >= tfxOvertime_velocity && type != tfxOvertime_velocity_adjuster;
	}

	bool IsGlobalGraph(GraphType type) {
		return type >= tfxGlobal_life && type <= tfxGlobal_splatter;
	}

	ParticleManager::~ParticleManager() {
	}

	EffectEmitter &ParticleManager::operator[] (unsigned int index) {
		return effects[current_ebuff][index];
	}
	
	void ParticleManager::AddEffect(EffectEmitter &effect, unsigned int buffer) {
		if (effects[buffer].current_size == effects[buffer].capacity)
			return;
		unsigned int parent_index = effects[buffer].current_size++;
		effects[buffer][parent_index] = effect;
		effects[buffer][parent_index].active_children = 0;
		effects[buffer][parent_index].flags &= ~tfxEmitterStateFlags_retain_matrix;
		effects[buffer][parent_index].pm = this;
		effects[buffer][parent_index].ResetParents();
		for (auto &e : effect.sub_effectors) {
			if (!FreeEffectCapacity())
				break;
			if (e.flags & tfxEmitterStateFlags_enabled) {
				unsigned int index = effects[buffer].current_size++;
				effects[buffer][index] = e;
				EffectEmitter &emitter = effects[buffer].back();
				emitter.parent = &effects[buffer][parent_index];
				emitter.next_ptr = emitter.parent;
				emitter.pm = this;
				emitter.active_children = 0;
				emitter.flags &= ~tfxEmitterStateFlags_retain_matrix;
				effects[buffer][parent_index].active_children++;
			}
		}
		effects[buffer][parent_index].NoTweenNextUpdate();
	}

	 void ParticleManager::AddEffect(EffectEmitterTemplate &effect, unsigned int buffer) {
		AddEffect(effect.effect_template, current_ebuff);
	}

	uint32_t ParticleManager::AddParticle(unsigned int layer, Particle &p) {
		assert(particles[layer][current_pbuff].current_size != particles[layer][current_pbuff].capacity);
		particles[layer][current_pbuff][particles[layer][current_pbuff].current_size] = p;
		particles[layer][current_pbuff].current_size++;
		return (uint32_t)particles[layer][current_pbuff].current_size - 1;
	}

	Particle& ParticleManager::GrabParticle(unsigned int layer) {
		assert(particles[layer][current_pbuff].current_size != particles[layer][current_pbuff].capacity);
		particles[layer][current_pbuff].current_size++;
		//particles[layer][current_pbuff][particles[layer][current_pbuff].current_size -1].particle_id = particle_id++;
		return particles[layer][current_pbuff][particles[layer][current_pbuff].current_size - 1];
	}

	void ParticleManager::Update() {
		unsigned int index = 0;

		unsigned int next_buffer = !current_ebuff;
		effects[next_buffer].clear();

		for (auto &e : effects[current_ebuff]) {
			e.Update();
			if (e.type == tfxEffect) {
				if (e.active_children > 0) {
					e.next_ptr = SetNextEffect(e, next_buffer);
				}
				else {
					e.next_ptr = nullptr;
					e.sub_effectors.free_all();
				}
			}
			else {
				if (e.timeout_counter < e.timeout) {
					e.next_ptr = SetNextEffect(e, next_buffer);
				}
				else {
					e.next_ptr = nullptr;
					e.sub_effectors.free_all();
				}
			}
			index++;
		}

		current_ebuff = next_buffer;

		next_buffer = !current_pbuff;

		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			particles[layer][next_buffer].clear();

			index = 0;
			for (auto &p : particles[layer][current_pbuff]) {
				p.parent = p.parent->next_ptr;

				if (!(p.flags & tfxParticleFlags_fresh)) {

					p.captured = p.world;

					if (ControlParticle(p, *p.parent)) {
						TransformParticle(p, *p.parent);

						p.next_ptr = SetNextParticle(layer, p, next_buffer);
					}
					else {
						p.parent->particle_count--;
						p.next_ptr = nullptr;
					}
				}
				else {
					p.flags &= ~tfxParticleFlags_fresh;
					p.next_ptr = SetNextParticle(layer, p, next_buffer);
				}
				index++;
			}
		}

		current_pbuff = next_buffer;

		update_base_values = false;

	}

	inline Particle* ParticleManager::SetNextParticle(unsigned int layer, Particle &p, unsigned int buffer) {
		unsigned int index = particles[layer][buffer].current_size++;
		assert(index < particles[layer][buffer].capacity);
		particles[layer][buffer][index] = p;
		return &particles[layer][buffer][index];
	}

	inline EffectEmitter* ParticleManager::SetNextEffect(EffectEmitter &e, unsigned int buffer) {
		unsigned int index = effects[buffer].current_size++;
		assert(index < effects[buffer].capacity);
		effects[buffer][index] = e;
		return &effects[buffer][index];
	}

	void ParticleManager::Render(float tween, void *data) {
		if (!render_func)
			return;
		for (auto &p : particles[current_pbuff]) {
			render_func(tween, &p, data);
		}
	}
	
	tfxvec<Particle> *ParticleManager::GetParticleBuffer(unsigned int layer) {
		return &particles[layer][current_pbuff];
	}
	
	tfxvec<EffectEmitter> *ParticleManager::GetEffectBuffer() {
		return &effects[current_ebuff];
	}

	void ParticleManager::SetRenderCallback(void func(float, void*, void*)) {
		render_func = func;
	}

	void ParticleManager::Init(unsigned int effects_limit, unsigned int particle_limit_per_layer) {
		max_particles_per_layer = particle_limit_per_layer;
		max_effects = effects_limit;

		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			particles[layer][0].resize(max_particles_per_layer);
			particles[layer][1].resize(max_particles_per_layer);
			particles[layer][0].clear();
			particles[layer][1].clear();
		}
		effects[0].create_pool(max_effects);
		effects[1].create_pool(max_effects);
		effects[0].clear();
		effects[1].clear();

	}

	uint32_t ParticleManager::ParticleCount() {
		unsigned int count = 0;
		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			count += particles[layer][current_pbuff].size();
		}
		return count;
	}

	void ParticleManager::ClearAll() {
		for (unsigned int layer = 0; layer != tfxLAYERS; ++layer) {
			particles[layer][0].clear();
			particles[layer][1].clear();
		}
		for (unsigned int i = 0; i != 2; ++i) {
			for (auto &e : effects[i]) {
				e.sub_effectors.free_all();
			}
			effects[i].clear();
		}
		particle_id = 0;
	}
	void ParticleManager::SoftExpireAll() {
		for (auto &e : effects[current_ebuff]) {
			e.flags |= tfxEmitterStateFlags_stop_spawning;
		}
	}
	void ParticleManager::ClearDepths() {
		for (unsigned int i = 0; i != 2; ++i) {
			for (auto &e : effects[i]) {
				e.sub_effectors.free_all();
			}
			effects[i].clear();
		}
	}
	void ParticleManager::SetLookUpMode(LookupMode mode) {
		if (mode == tfxPrecise) {
			lookup_overtime_callback = LookupPreciseOvertime;
			lookup_callback = LookupPrecise;
		}
		else {
			lookup_overtime_callback = LookupFastOvertime;
			lookup_callback = LookupFast;
		}
		lookup_mode = mode;
	}

	void ParticleManager::UpdateBaseValues() {
		update_base_values = true;
	}

	bool HasDataValue(tfxStorageMap<DataEntry> &config, tfxText key) {
		return config.ValidName(key);
	}

	void AddDataValue(tfxStorageMap<DataEntry> &map, tfxText key, const char *value) {
		DataEntry entry;
		entry.type = tfxString;
		entry.key = key;
		entry.str_value = value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<DataEntry> &map, tfxText key, int value) {
		DataEntry entry;
		entry.type = tfxSInt;
		entry.key = key;
		entry.int_value = value;
		entry.bool_value = (bool)value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<DataEntry> &map, tfxText key, bool value) {
		DataEntry entry;
		entry.type = tfxBool;
		entry.key = key;
		entry.bool_value = value;
		entry.int_value = (int)value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<DataEntry> &map, tfxText key, double value) {
		DataEntry entry;
		entry.type = tfxDouble;
		entry.key = key;
		entry.double_value = value;
		map.Insert(key, entry);
	}

	void AddDataValue(tfxStorageMap<DataEntry> &map, tfxText key, float value) {
		DataEntry entry;
		entry.type = tfxFloat;
		entry.key = key;
		entry.float_value = value;
		map.Insert(key, entry);
	}

	tfxText &GetDataStrValue(tfxStorageMap<DataEntry> &map, const char* key) {
		return map.At(key).str_value;
	}
	int& GetDataIntValue(tfxStorageMap<DataEntry> &map, const char* key) {
		return map.At(key).int_value;
	}
	float& GetDataFloatValue(tfxStorageMap<DataEntry> &map, const char* key) {
		return map.At(key).float_value;
	}

	void SaveDataFile(tfxStorageMap<DataEntry> &map, const char* path) {
		std::ofstream file(path);

		if (map.Size()) {
			for (auto &entry : map.data) {
				tfxText ini_line = entry.key;
				ini_line.Appendf("=");
				switch (entry.type) {
				case tfxString:
					ini_line.Appendf(entry.str_value.c_str());
					break;
				case tfxSInt:
					ini_line.Appendf("%i", entry.int_value);
					break;
				case tfxFloat:
					ini_line.Appendf("%f", entry.float_value);
					break;
				case tfxBool:
					ini_line.Appendf("%i", (int)entry.bool_value);
					break;
				}
				ini_line.Appendf("\n");
				file << ini_line.c_str();
			}
		}

		file.close();

	}

	void LoadDataFile(tfxStorageMap<DataEntry> &map, const char* path) {
		FILE* fp;
		fp = fopen(path, "r");
		if (fp == NULL) {
			return;
		}

		const size_t max_line_length = 256;
		char buffer[max_line_length];

		while (fgets(buffer, max_line_length, fp)) {
			buffer[strcspn(buffer, "\n")] = 0;
			tfxText str = buffer;
			tfxvec<tfxText> pair = SplitString(str, 61);
			if (pair.size() == 2) {
				tfxText key = pair[0];
				DataType t = data_types.eff.At(pair[0]);
				if (t == tfxBool) {
					AddDataValue(map, key, (bool)atoi(pair[1].c_str()));
				}
				if (t == tfxSInt) {
					AddDataValue(map, key, atoi(pair[1].c_str()));
				}
				else if(t == tfxFloat) {
					AddDataValue(map, key, (float)atof(pair[1].c_str()));
				}
				else if (t == tfxString) {
					AddDataValue(map, key, pair[1].c_str());
				}
			}
		}

		int close = fclose(fp);

	}

	tfxvec<tfxText> SplitString(const tfx::tfxText &str, char delim) {
		tfxvec<tfxText> ret;

		tfxText line;
		for (char c : str.string) {
			if (c == delim && line.Length() && c != NULL) {
				ret.push_back(line);
				line.Clear();
			}
			else if(c != NULL) {
				line.Append(c);
			}
		}

		if (line.Length()) {
			ret.push_back(line);
		}

		return ret;
	}

	bool StringIsUInt(const tfxText &s) {

		for (auto c : s.string) {
			if (!std::isdigit(c) && c != 0)
				return false;
		}

		return true;
	}

	int GetDataType(const tfxText &s) {
		if (s.Length() == 0)
			return tfxString;

		if (s[0] >= 48 && s[0] <= 57) {
			size_t size;
			auto int_check = std::stoi(s.c_str(), &size);
			if (size == s.Length())
				return tfxSInt;

			auto float_check = std::stof(s.c_str(), &size);
			if (size == s.Length())
				return tfxFloat;
		}

		return tfxString;
	}

	//Get a graph by GraphID
	Graph &GetGraph(EffectLibrary &library, GraphID &graph_id) {
		GraphType type = graph_id.type;

		if (type < tfxGlobalCount) {
			return ((Graph*)&library.global_graphs[graph_id.graph_id])[type];
		}
		else if (type >= tfxPropertyStart && type < tfxBaseStart) {
			int ref = type - tfxPropertyStart;
			return ((Graph*)&library.property_graphs[graph_id.graph_id])[ref];
		}
		else if (type >= tfxBaseStart && type < tfxVariationStart) {
			int ref = type - tfxBaseStart;
			return ((Graph*)&library.base_graphs[graph_id.graph_id])[ref];
		}
		else if (type >= tfxVariationStart && type < tfxOvertimeStart) {
			int ref = type - tfxVariationStart;
			return ((Graph*)&library.variation_graphs[graph_id.graph_id])[ref];
		}
		else if (type >= tfxOvertimeStart) {
			int ref = type - tfxOvertimeStart;
			return ((Graph*)&library.overtime_graphs[graph_id.graph_id])[ref];
		}

		assert(0);	//This function must return a value, make sure the graph_id is valid

		return((Graph*)&library.overtime_graphs[graph_id.graph_id])[type];

	}

	//Get a node by GraphID
	//Graph &GetGraphNode(EffectLibrary &library, GraphNodeID &graph_id) {
	//}


	//API Functions
	int GetShapesInLibrary(const char *filename) {
		tfxvec<EffectEmitter> effect_stack;
		int context = 0;
		int error = 0;
		int uid = 0;
		unsigned int current_global_graph = 0;

		mz_zip_archive zip_archive;
		memset(&zip_archive, 0, sizeof(zip_archive));

		bool status = mz_zip_reader_init_file(&zip_archive, filename, 0);
		if (!status)
		{
			error = -2;
		}

		if ((int)mz_zip_reader_get_num_files(&zip_archive) == 0)
			error = -3;

		size_t uncomp_size;
		std::stringstream d;
		void *data = mz_zip_reader_extract_file_to_heap(&zip_archive, "data.txt", &uncomp_size, 0);

		if (!data)
			error = -2;
		else
			d << (char*)data;


		if (error < 0) {
			mz_zip_reader_end(&zip_archive);
			return error;
		}

		int shape_count = 0;

		while (!d.eof()) {
			std::string line;
			std::getline(d, line);
			bool context_set = false;
			if (StringIsUInt(line.c_str()) && context != tfxStartShapes) {
				context = atoi(line.c_str());
				if (context == tfxEndShapes)
					break;
				context_set = true;
			}
			if (context_set == false) {
				tfxvec<tfxText> pair = SplitString(std::string(line).c_str());
				if (pair.size() != 2) {
					pair = SplitString(std::string(line).c_str(), 44);
					if (pair.size() < 2) {
						error = 1;
						break;
					}
				}
				if (context == tfxStartShapes) {
					if (pair.size() >= 5) {
						int frame_count = atoi(pair[2].c_str());
						shape_count += frame_count;
					}
				}
			}
		}

		mz_zip_reader_end(&zip_archive);
		return shape_count;
	}

	int LoadEffectLibrary(const char *filename, EffectLibrary &lib, void(*shape_loader)(const char* filename, ImageData &image_data, void *raw_image_data, int image_size, void *user_data), void *user_data) {
		assert(shape_loader);
		lib.Clear();

		tfxvec<EffectEmitter> effect_stack;
		int context = 0;
		int error = 0;
		int uid = 0;
		unsigned int current_global_graph = 0;

		mz_zip_archive zip_archive;
		memset(&zip_archive, 0, sizeof(zip_archive));

		bool status = mz_zip_reader_init_file(&zip_archive, filename, 0);
		if (!status)
		{
			error = -2;
		}

		if ((int)mz_zip_reader_get_num_files(&zip_archive) == 0)
			error = -3;

		size_t uncomp_size;
		std::stringstream d;
		void *data = mz_zip_reader_extract_file_to_heap(&zip_archive, "data.txt", &uncomp_size, 0);

		if (!data)
			error = -2;
		else 
			d << (char*)data;


		if (error < 0) {
			mz_free(data);
			mz_zip_reader_end(&zip_archive);
			return error;
		}

		while (!d.eof()) {
			std::string line;
			std::getline(d, line);
			bool context_set = false;

			if (StringIsUInt(line.c_str()) && context != tfxStartShapes) {
				context = atoi(line.c_str());
				if (context == tfxEndShapes)
					break;
				context_set = true;
				if (context == tfxStartEffect) {
					EffectEmitter effect;
					effect.library = &lib;
					if (effect_stack.size() == 0) { //Only root effects get the global graphs
						lib.AddEffectGraphs(effect);
						effect.ResetEffectGraphs(false);
						current_global_graph = effect.global;
					}
					effect.uid = uid++;
					effect.properties = EmitterProperties();
					effect.type = EffectEmitterType::tfxEffect;
					effect_stack.push_back(effect);
					lib.AddAnimationSettings(effect_stack.back());
				}
				if (context == tfxStartEmitter) {
					EffectEmitter emitter;
					emitter.library = &lib;
					lib.AddEmitterGraphs(emitter);
					emitter.uid = uid++;
					emitter.type = EffectEmitterType::tfxEmitter;
					emitter.ResetEmitterGraphs(false);
					effect_stack.push_back(emitter);
				}
			}

			if (context_set == false) {
				tfxvec<tfxText> pair = SplitString(line.c_str());
				if (pair.size() != 2) {
					pair = SplitString(line.c_str(), 44);
					if (pair.size() < 2) {
						error = 1;
						break;
					}
				}

				if (context == tfxStartEffect) {
					switch (data_types.eff.At(pair[0])) {
					case tfxUint:
						AssignEffectorProperty(effect_stack.back(), pair[0], (unsigned int)atoi(pair[1].c_str()));
						break;
					case tfxFloat:
						AssignEffectorProperty(effect_stack.back(), pair[0], (float)atof(pair[1].c_str()));
						break;
					case tfxSInt:
						AssignEffectorProperty(effect_stack.back(), pair[0], atoi(pair[1].c_str()));
						break;
					case tfxBool:
						AssignEffectorProperty(effect_stack.back(), pair[0], (bool)(atoi(pair[1].c_str())));
						break;
					case tfxString:
						AssignEffectorProperty(effect_stack.back(), pair[0], pair[1]);
						break;
					}
				}

				if (context == tfxStartAnimationSettings) {
					switch (data_types.eff.At(pair[0])) {
					case tfxUint:
						AssignEffectorProperty(effect_stack.back(), pair[0], (uint32_t)atoi(pair[1].c_str()));
						break;
					case tfxFloat:
						AssignEffectorProperty(effect_stack.back(), pair[0], (float)atof(pair[1].c_str()));
						break;
					case tfxSInt:
						AssignEffectorProperty(effect_stack.back(), pair[0], atoi(pair[1].c_str()));
						break;
					case tfxBool:
						AssignEffectorProperty(effect_stack.back(), pair[0], (bool)atoi(pair[1].c_str()));
						break;
					case tfxString:
						AssignEffectorProperty(effect_stack.back(), pair[0], pair[1]);
						break;
					}
				}

				if (context == tfxStartEmitter) {
					switch (data_types.eff.At(pair[0])) {
					case tfxUint:
						AssignEffectorProperty(effect_stack.back(), pair[0], (uint32_t)atoi(pair[1].c_str()));
						break;
					case tfxFloat:
						AssignEffectorProperty(effect_stack.back(), pair[0], (float)atof(pair[1].c_str()));
						break;
					case tfxSInt:
						AssignEffectorProperty(effect_stack.back(), pair[0], atoi(pair[1].c_str()));
						break;
					case tfxBool:
						AssignEffectorProperty(effect_stack.back(), pair[0], (bool)atoi(pair[1].c_str()));
						break;
					case tfxString:
						AssignEffectorProperty(effect_stack.back(), pair[0], pair[1]);
						break;
					}
				}

				if (context == tfxStartGraphs && effect_stack.back().type == tfxEmitter) {
					AssignGraphData(effect_stack.back(), pair);
				}
				else if (context == tfxStartGraphs && effect_stack.back().type == tfxEffect) {
					if (effect_stack.size() == 1)
						AssignGraphData(effect_stack.back(), pair);
				}

				if (context == tfxStartShapes) {
					if (pair.size() >= 5) {
						ShapeData s;
						strcpy_s(s.name, pair[0].c_str());
						s.shape_index = atoi(pair[1].c_str());
						s.frame_count = atoi(pair[2].c_str());
						s.width = atoi(pair[3].c_str());
						s.height = atoi(pair[4].c_str());
						if (pair.size() > 5)
							s.import_filter = atoi(pair[5].c_str());
						if (s.import_filter < 0 || s.import_filter>1)
							s.import_filter = 0;

						void *tmp = mz_zip_reader_extract_file_to_heap(&zip_archive, s.name, &uncomp_size, 0);
						ImageData image_data;
						image_data.animation_frames = (float)s.frame_count;
						image_data.image_size = tfxVec2((float)s.width, (float)s.height);
						image_data.shape_index = s.shape_index;
						image_data.import_filter = s.import_filter;

						if (tmp) {
							shape_loader(s.name, image_data, tmp, (int)uncomp_size, user_data);
						}

						if (!image_data.ptr) {
							uid = -5;
							break;
						}
						lib.particle_shapes.InsertByInt(s.shape_index, image_data);
						free(tmp);
					}
				}
			}

			if (context == tfxEndEmitter) {
				effect_stack.back().InitialiseUninitialisedGraphs();
				effect_stack.back().UpdateMaxLife();
				effect_stack.parent().sub_effectors.push_back(effect_stack.back());
				effect_stack.pop();
			}

			if (context == tfxEndEffect) {
				effect_stack.back().ReIndex();
				if (effect_stack.size() > 1)
					effect_stack.parent().sub_effectors.push_back(effect_stack.back());
				else {
					lib.effects.push_back(effect_stack.back());
					effect_stack.back().InitialiseUninitialisedGraphs();
				}
				effect_stack.pop();
			}

		}

		d.str(std::string());

		mz_free(data);
		mz_zip_reader_end(&zip_archive);

		//Returning anything over 0 means that effects were loaded ok
		//-1 = No effects were loaded
		//-2 = Failed to initialise the zip file
		//-3 = No files were found in the zip file
		//-4 = No data was found in the zip file
		//-5 = ShapeLoader did not add a pointer into the image data

		if (uid >= 0) {
			lib.CompileAllGraphs();
			lib.ReIndex();
			//lib.UpdateParticleShapeReferences(lib.effects, 1);
			lib.UpdateEffectPaths();
			lib.UpdateAllNodes();
			lib.SetMinMaxData();
		}

		return uid - 1;

	}

	 void StopSpawning(ParticleManager &pm) {
		pm.SoftExpireAll();
	}

	 void RemoveAllEffects(ParticleManager &pm) {
		pm.ClearAll();
	}

	 void InitParticleManager(ParticleManager &pm, unsigned int effects_limit, unsigned int particle_limit_per_layer) {
		pm.Init(effects_limit, particle_limit_per_layer);
	}

	void AddEffect(ParticleManager &pm, EffectEmitter &effect, float x, float y) {
		effect.Position(x, y);
		pm.AddEffect(effect, pm.current_ebuff);
	}

	void AddEffect(ParticleManager &pm, EffectEmitterTemplate &effect, float x, float y) {
		effect.effect_template.Position(x, y);
		pm.AddEffect(effect, pm.current_ebuff);
	}

	 void EffectEmitterTemplate::SetUserDataAll(void *data) {
		tfxvec<EffectEmitter*> stack;
		stack.push_back(&effect_template);
		while (stack.size()) {
			EffectEmitter *current = stack.pop_back();
			current->user_data = data;
			for (auto &sub : current->sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	 void EffectEmitterTemplate::SetUpdateCallbackAll(void(*update_callback)(EffectEmitter &effectemitter)) {
		tfxvec<EffectEmitter*> stack;
		stack.push_back(&effect_template);
		while (stack.size()) {
			EffectEmitter *current = stack.pop_back();
			current->update_callback = update_callback;
			for (auto &sub : current->sub_effectors) {
				stack.push_back(&sub);
			}
		}
	}

	void EffectEmitterTemplate::SetParticleUpdateCallback(tfxText path, void(*particle_update_callback)(Particle &particle)) {
		assert(paths.ValidName(path));
		EffectEmitter &e = *paths.At(path);
		assert(e.type == tfxEmitter);
		e.particle_update_callback = particle_update_callback;
	}

	void EffectEmitterTemplate::SetParticleOnSpawnCallback(tfxText path, void(*particle_onspawn_callback)(Particle &particle)) {
		assert(paths.ValidName(path));
		EffectEmitter &e = *paths.At(path);
		assert(e.type == tfxEmitter);
		e.particle_onspawn_callback = particle_onspawn_callback;
	}

}