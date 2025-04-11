from typing import List, Optional, Callable, NamedTuple
from BaseClasses import CollectionState
from .Options import TimespinnerOptions
from .PreCalculatedWeights import PreCalculatedWeights
from .LogicExtensions import TimespinnerLogic

EventId: Optional[int] = None


class LocationData(NamedTuple):
    region: str
    name: str
    code: Optional[int]
    rule: Optional[Callable[[CollectionState], bool]] = None


def get_location_datas(player: Optional[int], options: Optional[TimespinnerOptions],
                  precalculated_weights: Optional[PreCalculatedWeights]) -> List[LocationData]:
    flooded: Optional[PreCalculatedWeights] = precalculated_weights
    logic = TimespinnerLogic(player, options, precalculated_weights)

    # 1337000 - 1337155 Generic locations
    # 1337171 - 1337175 New Pickup checks
    # 1337246 - 1337249 Ancient Pyramid
    location_table: List[LocationData] = [
        # Present item locations
        LocationData('Tutorial', 'Tutorial: Yo Momma 1',  1337000),
        LocationData('Tutorial', 'Tutorial: Yo Momma 2',  1337001),
        LocationData('Lake desolation', 'Lake Desolation: Starter chest 2',  1337002),
        LocationData('Lake desolation', 'Lake Desolation: Starter chest 3',  1337003),
        LocationData('Lake desolation', 'Lake Desolation: Starter chest 1',  1337004),
        LocationData('Lake desolation', 'Lake Desolation (Lower): Timespinner Wheel room',  1337005),
        LocationData('Lake desolation', 'Lake Desolation: Forget me not chest',  1337006, lambda state: logic.has_fire(state) and state.can_reach('Upper Lake Serene', 'Region', player)),
        LocationData('Lake desolation', 'Lake Desolation (Lower): Chicken chest', 1337007, logic.has_timestop),
        LocationData('Lower lake desolation', 'Lake Desolation (Lower): Not so secret room',  1337008, logic.can_break_walls),
        LocationData('Lower lake desolation', 'Lake Desolation (Upper): Tank chest',  1337009, logic.has_timestop),
        LocationData('Upper lake desolation', 'Lake Desolation (Upper): Oxygen recovery room',  1337010),
        LocationData('Upper lake desolation', 'Lake Desolation (Upper): Secret room',  1337011, logic.can_break_walls),
        LocationData('Upper lake desolation', 'Lake Desolation (Upper): Double jump cave platform',  1337012, logic.has_doublejump),
        LocationData('Upper lake desolation', 'Lake Desolation (Upper): Double jump cave floor',  1337013),
        LocationData('Upper lake desolation', 'Lake Desolation (Upper): Sparrow chest',  1337014),
        LocationData('Upper lake desolation', 'Lake Desolation (Upper): Crash site pedestal',  1337015),
        LocationData('Upper lake desolation', 'Lake Desolation (Upper): Crash site chest 1',  1337016, lambda state: state.has('Killed Maw', player)),
        LocationData('Upper lake desolation', 'Lake Desolation (Upper): Crash site chest 2',  1337017, lambda state: state.has('Killed Maw', player)),
        LocationData('Eastern lake desolation', 'Lake Desolation: Kitty Boss',  1337018),
        LocationData('Library', 'Library: Basement',  1337019),
        LocationData('Library', 'Library: Warp gate',  1337020),
        LocationData('Library', 'Library: Librarian',  1337021),
        LocationData('Library', 'Library: Reading nook chest',  1337022),
        LocationData('Library', 'Library: Storage room chest 1',  1337023, logic.has_keycard_D),
        LocationData('Library', 'Library: Storage room chest 2',  1337024, logic.has_keycard_D),
        LocationData('Library', 'Library: Storage room chest 3',  1337025, logic.has_keycard_D),
        LocationData('Library top', 'Library: Backer room chest 5',  1337026),
        LocationData('Library top', 'Library: Backer room chest 4',  1337027),
        LocationData('Library top', 'Library: Backer room chest 3',  1337028),
        LocationData('Library top', 'Library: Backer room chest 2',  1337029),
        LocationData('Library top', 'Library: Backer room chest 1',  1337030),
        LocationData('Varndagroth tower left', 'Varndagroth Towers (Left): Elevator Key not required',  1337031),
        LocationData('Varndagroth tower left', 'Varndagroth Towers (Left): Ye olde Timespinner',  1337032),
        LocationData('Varndagroth tower left', 'Varndagroth Towers (Left): Bottom floor',  1337033, logic.has_keycard_C),
        LocationData('Varndagroth tower left', 'Varndagroth Towers (Left): Air vents secret',  1337034, logic.can_break_walls),
        LocationData('Varndagroth tower left', 'Varndagroth Towers (Left): Elevator chest',  1337035, lambda state: state.has('Elevator Keycard', player)),
        LocationData('Varndagroth tower right (upper)', 'Varndagroth Towers: Bridge',  1337036),
        LocationData('Varndagroth tower right (elevator)', 'Varndagroth Towers (Right): Elevator chest',  1337037),
        LocationData('Varndagroth tower right (upper)', 'Varndagroth Towers (Right): Elevator card chest',  1337038, lambda state: state.has('Elevator Keycard', player) or logic.has_doublejump(state)),
        LocationData('Varndagroth tower right (upper)', 'Varndagroth Towers (Right): Air vents right chest',  1337039, lambda state: state.has('Elevator Keycard', player) or logic.has_doublejump(state)),
        LocationData('Varndagroth tower right (upper)', 'Varndagroth Towers (Right): Air vents left chest',  1337040, lambda state: state.has('Elevator Keycard', player) or logic.has_doublejump(state)),
        LocationData('Varndagroth tower right (lower)', 'Varndagroth Towers (Right): Bottom floor',  1337041),
        LocationData('Varndagroth tower right (elevator)', 'Varndagroth Towers (Right): Varndagroth',  1337042, logic.has_keycard_C),
        LocationData('Varndagroth tower right (elevator)', 'Varndagroth Towers (Right): Spider Hell',  1337043, logic.has_keycard_A),
        LocationData('Skeleton Shaft', 'Sealed Caves (Xarion): Skeleton',  1337044),
        LocationData('Sealed Caves (Xarion)', 'Sealed Caves (Xarion): Shroom jump room',  1337045, logic.has_timestop),
        LocationData('Sealed Caves (Xarion)', 'Sealed Caves (Xarion): Double shroom room',  1337046),
        LocationData('Sealed Caves (Xarion)', 'Sealed Caves (Xarion): Jacksquat room',  1337047, logic.has_forwarddash_doublejump),
        LocationData('Sealed Caves (Xarion)', 'Sealed Caves (Xarion): Below Jacksquat room',  1337048),
        LocationData('Sealed Caves (Xarion)', 'Sealed Caves (Xarion): Secret room',  1337049, logic.can_break_walls),
        LocationData('Sealed Caves (Xarion)', 'Sealed Caves (Xarion): Bottom left room',  1337050),
        LocationData('Sealed Caves (Xarion)', 'Sealed Caves (Xarion): Last chance before Xarion',  1337051, logic.has_doublejump),
        LocationData('Sealed Caves (Xarion)', 'Sealed Caves (Xarion): Xarion',  1337052, lambda state: not flooded.flood_xarion or state.has('Water Mask', player)),
        LocationData('Sealed Caves (Sirens)', 'Sealed Caves (Sirens): Water hook',  1337053, lambda state: state.has('Water Mask', player)),
        LocationData('Sealed Caves (Sirens)', 'Sealed Caves (Sirens): Siren room underwater right',  1337054, lambda state: state.has('Water Mask', player)),
        LocationData('Sealed Caves (Sirens)', 'Sealed Caves (Sirens): Siren room underwater left',  1337055, lambda state: state.has('Water Mask', player)),
        LocationData('Sealed Caves (Sirens)', 'Sealed Caves (Sirens): Cave after sirens chest 1',  1337056),
        LocationData('Sealed Caves (Sirens)', 'Sealed Caves (Sirens): Cave after sirens chest 2',  1337057),
        LocationData('Military Fortress', 'Military Fortress: Bomber chest',  1337058, lambda state: state.has('Timespinner Wheel', player) and logic.has_doublejump_of_npc(state)),
        LocationData('Military Fortress', 'Military Fortress: Close combat room',  1337059),
        LocationData('Military Fortress (hangar)', 'Military Fortress: Soldiers bridge',  1337060),
        LocationData('Military Fortress (hangar)', 'Military Fortress: Giantess room',  1337061),
        LocationData('Military Fortress (hangar)', 'Military Fortress: Giantess bridge',  1337062),
        LocationData('Military Fortress (hangar)', 'Military Fortress: B door chest 2',  1337063, lambda state: logic.has_keycard_B(state) and (state.has('Water Mask', player) if flooded.flood_lab else logic.has_doublejump(state))),
        LocationData('Military Fortress (hangar)', 'Military Fortress: B door chest 1',  1337064, lambda state: logic.has_keycard_B(state) and (state.has('Water Mask', player) if flooded.flood_lab else logic.has_doublejump(state))),
        LocationData('Military Fortress (hangar)', 'Military Fortress: Pedestal',  1337065, lambda state: state.has('Water Mask', player) if flooded.flood_lab else (logic.has_doublejump_of_npc(state) or logic.has_forwarddash_doublejump(state))),
        LocationData('The lab', 'Lab: Coffee break',  1337066),
        LocationData('The lab', 'Lab: Lower trash right',  1337067, logic.has_doublejump),
        LocationData('The lab', 'Lab: Lower trash left',  1337068, lambda state: logic.has_doublejump_of_npc(state) if options.lock_key_amadeus else logic.has_upwarddash ),
        LocationData('The lab', 'Lab: Below lab entrance',  1337069, logic.has_doublejump),
        LocationData('The lab (power off)', 'Lab: Trash jump room',  1337070, lambda state: not options.lock_key_amadeus or logic.has_doublejump_of_npc(state) ),
        LocationData('The lab (power off)', 'Lab: Dynamo Works',  1337071, lambda state: not options.lock_key_amadeus or (state.has_all(('Lab Access Research', 'Lab Access Dynamo'), player)) ),
        LocationData('The lab (upper)', 'Lab: Genza (Blob Mom)',  1337072),
        LocationData('The lab (power off)', 'Lab: Experiment #13',  1337073, lambda state: not options.lock_key_amadeus or state.has('Lab Access Experiment', player) ),
        LocationData('The lab (upper)', 'Lab: Download and chest room chest',  1337074),
        LocationData('The lab (upper)', 'Lab: Lab secret',  1337075, logic.can_break_walls),
        LocationData('The lab (power off)', 'Lab: Spider Hell',  1337076, lambda state: logic.has_keycard_A and not options.lock_key_amadeus or state.has('Lab Access Research', player)),
        LocationData('Emperors tower', 'Emperor\'s Tower: Courtyard bottom chest',  1337077),
        LocationData('Emperors tower', 'Emperor\'s Tower: Courtyard floor secret',  1337078, lambda state: logic.has_upwarddash(state) and logic.can_break_walls(state)),
        LocationData('Emperors tower', 'Emperor\'s Tower: Courtyard upper chest',  1337079, lambda state: logic.has_upwarddash(state)),
        LocationData('Emperors tower', 'Emperor\'s Tower: Galactic sage room',  1337080),
        LocationData('Emperors tower', 'Emperor\'s Tower: Bottom right tower',  1337081),
        LocationData('Emperors tower', 'Emperor\'s Tower: Wayyyy up there',  1337082, logic.has_doublejump_of_npc),
        LocationData('Emperors tower', 'Emperor\'s Tower: Left tower balcony',  1337083),
        LocationData('Emperors tower', 'Emperor\'s Tower: Emperor\'s Chambers chest',  1337084),
        LocationData('Emperors tower', 'Emperor\'s Tower: Emperor\'s Chambers pedestal',  1337085),
        LocationData('Emperors tower', 'Killed Emperor', EventId),

        # Past item locations
        LocationData('Refugee Camp', 'Refugee Camp: Neliste\'s Bra',  1337086),
        LocationData('Refugee Camp', 'Refugee Camp: Storage chest 3',  1337087),
        LocationData('Refugee Camp', 'Refugee Camp: Storage chest 2',  1337088),
        LocationData('Refugee Camp', 'Refugee Camp: Storage chest 1',  1337089),
        LocationData('Forest', 'Forest: Refugee camp roof',  1337090),
        LocationData('Forest', 'Forest: Bat jump ledge',  1337091, lambda state: logic.has_doublejump_of_npc(state) or logic.has_forwarddash_doublejump(state) or logic.has_fastjump_on_npc(state)),
        LocationData('Forest', 'Forest: Green platform secret',  1337092, logic.can_break_walls),
        LocationData('Forest', 'Forest: Rats guarded chest',  1337093),
        LocationData('Forest', 'Forest: Waterfall chest 1',  1337094, lambda state: state.has('Water Mask', player)),
        LocationData('Forest', 'Forest: Waterfall chest 2',  1337095, lambda state: state.has('Water Mask', player)),
        LocationData('Forest', 'Forest: Batcave',  1337096),
        LocationData('Forest', 'Castle Ramparts: In the moat',  1337097, lambda state: not flooded.flood_moat or state.has('Water Mask', player)),
        LocationData('Left Side forest Caves', 'Forest: Before Serene single bat cave',  1337098),
        LocationData('Upper Lake Serene', 'Lake Serene (Upper): Rat nest',  1337099),
        LocationData('Upper Lake Serene', 'Lake Serene (Upper): Double jump cave platform',  1337100, logic.has_doublejump),
        LocationData('Upper Lake Serene', 'Lake Serene (Upper): Double jump cave floor',  1337101),
        LocationData('Upper Lake Serene', 'Lake Serene (Upper): Cave secret',  1337102, logic.can_break_walls),
        LocationData('Upper Lake Serene', 'Lake Serene: Before Big Bird', 1337175),
        LocationData('Upper Lake Serene', 'Lake Serene: Behind the vines',  1337103),
        LocationData('Upper Lake Serene', 'Lake Serene: Pyramid keys room',  1337104),
        LocationData('Upper Lake Serene', 'Lake Serene (Upper): Chicken ledge', 1337174),
        LocationData('Lower Lake Serene', 'Lake Serene (Lower): Deep dive',  1337105),
        LocationData('Left Side forest Caves', 'Lake Serene (Lower): Under the eels',  1337106, lambda state: state.has('Water Mask', player)),
        LocationData('Lower Lake Serene', 'Lake Serene (Lower): Water spikes room',  1337107),
        LocationData('Lower Lake Serene', 'Lake Serene (Lower): Underwater secret',  1337108, logic.can_break_walls),
        LocationData('Lower Lake Serene', 'Lake Serene (Lower): T chest',  1337109, lambda state: flooded.flood_lake_serene or logic.has_doublejump_of_npc(state)),
        LocationData('Left Side forest Caves', 'Lake Serene (Lower): Past the eels',  1337110, lambda state: state.has('Water Mask', player)),
        LocationData('Lower Lake Serene', 'Lake Serene (Lower): Underwater pedestal',  1337111, lambda state: flooded.flood_lake_serene or logic.has_doublejump(state)),
        LocationData('Caves of Banishment (upper)', 'Caves of Banishment (Maw): Shroom jump room',  1337112, lambda state: flooded.flood_maw or logic.has_doublejump(state)),
        LocationData('Caves of Banishment (upper)', 'Caves of Banishment (Maw): Secret room',  1337113, lambda state: logic.can_break_walls(state) and (not flooded.flood_maw or state.has('Water Mask', player))),
        LocationData('Caves of Banishment (upper)', 'Caves of Banishment (Maw): Bottom left room',  1337114, lambda state: not flooded.flood_maw or state.has('Water Mask', player)),
        LocationData('Caves of Banishment (upper)', 'Caves of Banishment (Maw): Single shroom room',  1337115),
        LocationData('Caves of Banishment (upper)', 'Caves of Banishment (Maw): Jackpot room chest 1',  1337116, lambda state: flooded.flood_maw or logic.has_forwarddash_doublejump(state)),
        LocationData('Caves of Banishment (upper)', 'Caves of Banishment (Maw): Jackpot room chest 2',  1337117, lambda state: flooded.flood_maw or logic.has_forwarddash_doublejump(state)),
        LocationData('Caves of Banishment (upper)', 'Caves of Banishment (Maw): Jackpot room chest 3',  1337118, lambda state: flooded.flood_maw or logic.has_forwarddash_doublejump(state)),
        LocationData('Caves of Banishment (upper)', 'Caves of Banishment (Maw): Jackpot room chest 4',  1337119, lambda state: flooded.flood_maw or logic.has_forwarddash_doublejump(state)),
        LocationData('Caves of Banishment (upper)', 'Caves of Banishment (Maw): Pedestal',  1337120, lambda state: not flooded.flood_maw or state.has('Water Mask', player)),
        LocationData('Caves of Banishment (Maw)', 'Caves of Banishment (Maw): Last chance before Maw',  1337121, lambda state: state.has('Water Mask', player) if flooded.flood_maw else logic.has_doublejump(state)),
        LocationData('Caves of Banishment (Maw)', 'Caves of Banishment (Maw): Plasma Crystal', 1337173, lambda state: state.has_any({'Gas Mask', 'Talaria Attachment'}, player) and (not flooded.flood_maw or state.has('Water Mask', player))),
        LocationData('Caves of Banishment (Maw)', 'Killed Maw',  EventId, lambda state: state.has('Gas Mask', player) and (not flooded.flood_maw or state.has('Water Mask', player))),
        LocationData('Caves of Banishment (Maw)', 'Caves of Banishment (Maw): Mineshaft',  1337122, lambda state: state.has_any({'Gas Mask', 'Talaria Attachment'}, player) and (not flooded.flood_maw or state.has('Water Mask', player))),
        LocationData('Caves of Banishment (Sirens)', 'Caves of Banishment (Sirens): Wyvern room',  1337123),
        LocationData('Caves of Banishment (Sirens)', 'Caves of Banishment (Sirens): Siren room above water chest',  1337124),
        LocationData('Caves of Banishment (Sirens)', 'Caves of Banishment (Sirens): Siren room underwater left chest',  1337125, lambda state: state.has('Water Mask', player)),
        LocationData('Caves of Banishment (Sirens)', 'Caves of Banishment (Sirens): Siren room underwater right chest',  1337126, lambda state: state.has('Water Mask', player)),
        LocationData('Caves of Banishment (Sirens)', 'Caves of Banishment (Sirens): Siren room underwater right ground', 1337172, lambda state: state.has('Water Mask', player)),
        LocationData('Caves of Banishment (Sirens)', 'Caves of Banishment (Sirens): Water hook',  1337127, lambda state: state.has('Water Mask', player)),
        LocationData('Castle Ramparts', 'Castle Ramparts: Bomber chest',  1337128, logic.has_multiple_small_jumps_of_npc),
        LocationData('Castle Ramparts', 'Castle Ramparts: Freeze the engineer',  1337129, lambda state: state.has('Talaria Attachment', player) or logic.has_timestop(state)),
        LocationData('Castle Ramparts', 'Castle Ramparts: Giantess guarded room',  1337130),
        LocationData('Castle Ramparts', 'Castle Ramparts: Knight and archer guarded room',  1337131),
        LocationData('Castle Ramparts', 'Castle Ramparts: Pedestal',  1337132),
        LocationData('Castle Basement', 'Castle Basement: Secret pedestal',  1337133, logic.can_break_walls),
        LocationData('Castle Basement', 'Castle Basement: Clean the castle basement',  1337134),
        LocationData('Royal towers (lower)', 'Castle Keep: Yas queen room',  1337135, logic.has_pink),
        LocationData('Castle Basement', 'Castle Basement: Giantess guarded chest',  1337136),
        LocationData('Castle Basement', 'Castle Basement: Omelette chest',  1337137),
        LocationData('Castle Basement', 'Castle Basement: Just an egg',  1337138),
        LocationData('Castle Keep', 'Castle Keep: Under the twins',  1337139),
        LocationData('Castle Keep', 'Killed Twins',  EventId, logic.has_timestop),
        LocationData('Castle Keep', 'Castle Keep: Advisor jump', 1337171, logic.has_timestop),
        LocationData('Castle Keep', 'Castle Keep: Twins',  1337140, logic.has_timestop),
        LocationData('Castle Keep', 'Castle Keep: Royal guard tiny room',  1337141, lambda state: logic.has_doublejump(state) or logic.has_fastjump_on_npc(state)),
        LocationData('Royal towers (lower)', 'Royal Towers: Floor secret',  1337142, lambda state: logic.has_doublejump(state) and logic.can_break_walls(state)),
        LocationData('Royal towers', 'Royal Towers: Pre-climb gap',  1337143),
        LocationData('Royal towers', 'Royal Towers: Long balcony',  1337144, lambda state: not flooded.flood_courtyard or state.has('Water Mask', player)),
        LocationData('Royal towers', 'Royal Towers: Past bottom struggle juggle',  1337145, lambda state: flooded.flood_courtyard or logic.has_doublejump_of_npc(state)),
        LocationData('Royal towers', 'Royal Towers: Bottom struggle juggle',  1337146, logic.has_doublejump_of_npc),
        LocationData('Royal towers (upper)', 'Royal Towers: Top struggle juggle',  1337147, logic.has_doublejump_of_npc),
        LocationData('Royal towers (upper)', 'Royal Towers: No struggle required',  1337148),
        LocationData('Royal towers', 'Royal Towers: Right tower freebie',  1337149),
        LocationData('Royal towers (upper)', 'Royal Towers: Left tower small balcony',  1337150),
        LocationData('Royal towers (upper)', 'Royal Towers: Left tower royal guard',  1337151),
        LocationData('Royal towers (upper)', 'Royal Towers: Before Aelana',  1337152),
        LocationData('Royal towers (upper)', 'Killed Aelana',  EventId),
        LocationData('Royal towers (upper)', 'Royal Towers: Aelana\'s attic',  1337153, logic.has_upwarddash),
        LocationData('Royal towers (upper)', 'Royal Towers: Aelana\'s chest',  1337154),
        LocationData('Royal towers (upper)', 'Royal Towers: Aelana\'s pedestal',  1337155),

        # Ancient pyramid locations
        LocationData('Ancient Pyramid (entrance)', 'Ancient Pyramid: Why not it\'s right there',  1337246),
        LocationData('Ancient Pyramid (left)', 'Ancient Pyramid: Conviction guarded room',  1337247),
        LocationData('Ancient Pyramid (left)', 'Ancient Pyramid: Pit secret room',  1337248, lambda state: logic.can_break_walls(state) and (not flooded.flood_pyramid_shaft or state.has('Water Mask', player))),
        LocationData('Ancient Pyramid (left)', 'Ancient Pyramid: Regret chest',  1337249, lambda state: logic.can_break_walls(state) and (state.has('Water Mask', player) if flooded.flood_pyramid_shaft else logic.has_doublejump(state))),
        LocationData('Ancient Pyramid (right)', 'Ancient Pyramid: Nightmare Door chest',  1337236, lambda state: not flooded.flood_pyramid_back or state.has('Water Mask', player)),
        LocationData('Ancient Pyramid (right)', 'Killed Nightmare', EventId, lambda state: state.has_all({'Timespinner Wheel', 'Timespinner Spindle', 'Timespinner Gear 1', 'Timespinner Gear 2', 'Timespinner Gear 3'}, player) and (not flooded.flood_pyramid_back or state.has('Water Mask', player)))
    ]

    # 1337156 - 1337170 Downloads
    if not options or options.downloadable_items:
        location_table += ( 
            LocationData('Library', 'Library: Terminal 2 (Lachiem)',  1337156, lambda state: state.has('Tablet', player)),
            LocationData('Library', 'Library: Terminal 1 (Windaria)',  1337157, lambda state: state.has('Tablet', player)),
            # 1337158 Is lost in time
            LocationData('Library', 'Library: Terminal 3 (Emperor Nuvius)',  1337159, lambda state: state.has('Tablet', player)),
            LocationData('Library', 'Library: V terminal 1 (War of the Sisters)',  1337160, lambda state: state.has_all({'Tablet', 'Library Keycard V'}, player)),
            LocationData('Library', 'Library: V terminal 2 (Lake Desolation Map)',  1337161, lambda state: state.has_all({'Tablet', 'Library Keycard V'}, player)),
            LocationData('Library', 'Library: V terminal 3 (Vilete)',  1337162, lambda state: state.has_all({'Tablet', 'Library Keycard V'}, player)),
            LocationData('Library top', 'Library: Backer room terminal (Vandagray Metropolis Map)',  1337163, lambda state: state.has('Tablet', player)),
            LocationData('Varndagroth tower right (elevator)', 'Varndagroth Towers (Right): Medbay terminal (Bleakness Research)',  1337164, lambda state: state.has('Tablet', player) and logic.has_keycard_B(state)),
            LocationData('The lab (upper)', 'Lab: Download and chest room terminal (Experiment #13)',  1337165, lambda state: state.has('Tablet', player)),
            LocationData('The lab (power off)', 'Lab: Middle terminal (Amadeus Laboratory Map)',  1337166, lambda state: state.has('Tablet', player) and (not options.lock_key_amadeus or state.has('Lab Access Research', player))),
            LocationData('The lab (power off)', 'Lab: Sentry platform terminal (Origins)',  1337167, lambda state: state.has('Tablet', player) and (not options.lock_key_amadeus or state.has('Lab Access Genza', player) or logic.can_teleport_to(state, "Time", "GateDadsTower"))),
            LocationData('The lab', 'Lab: Experiment 13 terminal (W.R.E.C Farewell)',  1337168, lambda state: state.has('Tablet', player)),
            LocationData('The lab', 'Lab: Left terminal (Biotechnology)',  1337169, lambda state: state.has('Tablet', player)),
            LocationData('The lab (power off)', 'Lab: Right terminal (Experiment #11)',  1337170, lambda state: state.has('Tablet', player) and (not options.lock_key_amadeus or state.has('Lab Access Research', player)))
        )

    # 1337176 - 1337176 Cantoran
    if not options or options.cantoran:
        location_table += (
            LocationData('Left Side forest Caves', 'Lake Serene: Cantoran',  1337176),
        )

    # 1337177 - 1337198 Lore Checks
    if not options or options.lore_checks:
        location_table += (
            LocationData('Lower lake desolation', 'Lake Desolation: Memory - Coyote Jump (Time Messenger)',  1337177),
            LocationData('Library', 'Library: Memory - Waterway (A Message)',  1337178),
            LocationData('Library top', 'Library: Memory - Library Gap (Lachiemi Sun)',  1337179),
            LocationData('Library top', 'Library: Memory - Mr. Hat Portrait (Moonlit Night)',  1337180),
            LocationData('Varndagroth tower left', 'Varndagroth Towers (Left): Memory - Elevator (Nomads)',  1337181, lambda state: state.has('Elevator Keycard', player)),
            LocationData('Varndagroth tower right (lower)', 'Varndagroth Towers: Memory - Siren Elevator (Childhood)',  1337182, logic.has_keycard_B),
            LocationData('Varndagroth tower right (lower)', 'Varndagroth Towers (Right): Memory - Bottom (Faron)',  1337183),
            LocationData('Military Fortress', 'Military Fortress: Memory - Bomber Climb (A Solution)',  1337184, lambda state: state.has('Timespinner Wheel', player) and logic.has_doublejump_of_npc(state)),
            LocationData('The lab', 'Lab: Memory - Genza\'s Secret Stash 1 (An Old Friend)',  1337185, logic.can_break_walls),
            LocationData('The lab', 'Lab: Memory - Genza\'s Secret Stash 2 (Twilight Dinner)',  1337186, logic.can_break_walls),
            LocationData('Emperors tower', 'Emperor\'s Tower: Memory - Way Up There (Final Circle)',  1337187, logic.has_doublejump_of_npc),
            LocationData('Forest', 'Forest: Journal - Rats (Lachiem Expedition)',  1337188),
            LocationData('Forest', 'Forest: Journal - Bat Jump Ledge (Peace Treaty)',  1337189, lambda state: logic.has_doublejump_of_npc(state) or logic.has_forwarddash_doublejump(state) or logic.has_fastjump_on_npc(state)),
            LocationData('Forest', 'Forest: Journal - Floating in Moat (Prime Edicts)',  1337190, lambda state: not flooded.flood_moat or state.has('Water Mask', player)),
            LocationData('Castle Ramparts', 'Castle Ramparts: Journal - Archer + Knight (Declaration of Independence)',  1337191),
            LocationData('Castle Keep', 'Castle Keep: Journal - Under the Twins (Letter of Reference)',  1337192),
            LocationData('Castle Basement', 'Castle Basement: Journal - Castle Loop Giantess (Political Advice)',  1337193),
            LocationData('Royal towers (lower)', 'Royal Towers: Journal - Aelana\'s Room (Diplomatic Missive)',  1337194, logic.has_pink),
            LocationData('Royal towers (upper)', 'Royal Towers: Journal - Top Struggle Juggle Base (War of the Sisters)',  1337195),
            LocationData('Royal towers (upper)', 'Royal Towers: Journal - Aelana Boss (Stained Letter)',  1337196),
            LocationData('Royal towers', 'Royal Towers: Journal - Near Bottom Struggle Juggle (Mission Findings)',  1337197, lambda state: flooded.flood_courtyard or logic.has_doublejump_of_npc(state)),
            LocationData('Caves of Banishment (Maw)', 'Caves of Banishment (Maw): Journal - Lower Left Caves (Naivety)',  1337198, lambda state: not flooded.flood_maw or state.has('Water Mask', player))
        )

    # 1337199 - 1337232 Reserved for future use

    # 1337233 - 1337235 Pyramid Start checks
    if not options or options.pyramid_start: 
        location_table += (
            LocationData('Ancient Pyramid (entrance)', 'Dark Forest: Training Dummy',  1337233),
            LocationData('Ancient Pyramid (entrance)', 'Temporal Gyre: Forest Entrance',  1337234, lambda state: logic.has_upwarddash(state) or logic.can_teleport_to(state, "Time", "GateGyre")),
            LocationData('Ancient Pyramid (entrance)', 'Ancient Pyramid: Rubble',  1337235),
        )

    # 1337236 Nightmare door

    # 1337237 - 1337245 GyreArchives
    if not options or options.gyre_archives:
        location_table += (
            LocationData('Ravenlord\'s Lair', 'Ravenlord: Post fight (pedestal)',  1337237),
            LocationData('Ifrit\'s Lair', 'Ifrit: Post fight (pedestal)',  1337238),
            LocationData('Temporal Gyre', 'Temporal Gyre: Chest 1',  1337239),
            LocationData('Temporal Gyre', 'Temporal Gyre: Chest 2',  1337240),
            LocationData('Temporal Gyre', 'Temporal Gyre: Chest 3',  1337241),
            LocationData('Ravenlord\'s Lair', 'Ravenlord: Pre fight',  1337242),
            LocationData('Ravenlord\'s Lair', 'Ravenlord: Post fight (chest)',  1337243),
            LocationData('Ifrit\'s Lair', 'Ifrit: Pre fight',  1337244),
            LocationData('Ifrit\'s Lair', 'Ifrit: Post fight (chest)', 1337245),
        )

    # 1337250-1337781 Lanterns
    if not options or options.gyre_archives:
        location_table += (
            LocationData('<areaName>', 'ItemKey(1, 11, 106, 221)', 1337250)
            LocationData('<areaName>', 'ItemKey(1, 16, 106, 461)', 1337251)
            LocationData('<areaName>', 'ItemKey(1, 16, 282, 269)', 1337252)
            LocationData('<areaName>', 'ItemKey(1, 17, 282, 141)', 1337253)
            LocationData('<areaName>', 'ItemKey(1, 19, 218, 365)', 1337254)
            LocationData('<areaName>', 'ItemKey(1, 20, 250, 237)', 1337255)
            LocationData('<areaName>', 'ItemKey(1, 6, 474, 413)', 1337256)
            LocationData('<areaName>', 'ItemKey(1, 6, 570, 173)', 1337257)
            LocationData('<areaName>', 'ItemKey(1, 9, 122, 189)', 1337258)
            LocationData('<areaName>', 'ItemKey(1, 9, 314, 141)', 1337259)
            LocationData('<areaName>', 'ItemKey(10, 1, 136, 156)', 1337260)
            LocationData('<areaName>', 'ItemKey(10, 1, 264, 156)', 1337261)
            LocationData('<areaName>', 'ItemKey(10, 18, 120, 156)', 1337262)
            LocationData('<areaName>', 'ItemKey(10, 3, 232, 508)', 1337263)
            LocationData('<areaName>', 'ItemKey(10, 3, 232, 876)', 1337264)
            LocationData('<areaName>', 'ItemKey(10, 3, 568, 508)', 1337265)
            LocationData('<areaName>', 'ItemKey(10, 3, 568, 876)', 1337266)
            LocationData('<areaName>', 'ItemKey(10, 4, 1512, 140)', 1337267)
            LocationData('<areaName>', 'ItemKey(10, 4, 264, 140)', 1337268)
            LocationData('<areaName>', 'ItemKey(10, 4, 744, 140)', 1337269)
            LocationData('<areaName>', 'ItemKey(10, 5, 120, 556)', 1337270)
            LocationData('<areaName>', 'ItemKey(10, 5, 280, 556)', 1337271)
            LocationData('<areaName>', 'ItemKey(10, 6, 1048, 556)', 1337272)
            LocationData('<areaName>', 'ItemKey(10, 6, 152, 556)', 1337273)
            LocationData('<areaName>', 'ItemKey(10, 8, 1512, 140)', 1337274)
            LocationData('<areaName>', 'ItemKey(10, 8, 248, 140)', 1337275)
            LocationData('<areaName>', 'ItemKey(10, 8, 632, 140)', 1337276)
            LocationData('<areaName>', 'ItemKey(11, 0, 136, 169)', 1337277)
            LocationData('<areaName>', 'ItemKey(11, 0, 264, 169)', 1337278)
            LocationData('<areaName>', 'ItemKey(11, 16, 504, 169)', 1337279)
            LocationData('<areaName>', 'ItemKey(11, 16, 696, 169)', 1337280)
            LocationData('<areaName>', 'ItemKey(11, 17, 216, 281)', 1337281)
            LocationData('<areaName>', 'ItemKey(11, 18, 126, 144)', 1337282)
            LocationData('<areaName>', 'ItemKey(11, 18, 782, 112)', 1337283)
            LocationData('<areaName>', 'ItemKey(11, 19, 120, 169)', 1337284)
            LocationData('<areaName>', 'ItemKey(11, 19, 264, 169)', 1337285)
            LocationData('<areaName>', 'ItemKey(11, 22, 104, 169)', 1337286)
            LocationData('<areaName>', 'ItemKey(11, 22, 104, 505)', 1337287)
            LocationData('<areaName>', 'ItemKey(11, 22, 280, 169)', 1337288)
            LocationData('<areaName>', 'ItemKey(11, 22, 280, 505)', 1337289)
            LocationData('<areaName>', 'ItemKey(11, 23, 104, 169)', 1337290)
            LocationData('<areaName>', 'ItemKey(11, 23, 280, 169)', 1337291)
            LocationData('<areaName>', 'ItemKey(11, 3, 1448, 489)', 1337292)
            LocationData('<areaName>', 'ItemKey(11, 3, 152, 489)', 1337293)
            LocationData('<areaName>', 'ItemKey(11, 30, 152, 169)', 1337294)
            LocationData('<areaName>', 'ItemKey(11, 34, 136, 169)', 1337295)
            LocationData('<areaName>', 'ItemKey(11, 34, 264, 169)', 1337296)
            LocationData('<areaName>', 'ItemKey(11, 35, 120, 1209)', 1337297)
            LocationData('<areaName>', 'ItemKey(11, 35, 120, 169)', 1337298)
            LocationData('<areaName>', 'ItemKey(11, 35, 120, 537)', 1337299)
            LocationData('<areaName>', 'ItemKey(11, 35, 120, 905)', 1337300)
            LocationData('<areaName>', 'ItemKey(11, 35, 264, 169)', 1337301)
            LocationData('<areaName>', 'ItemKey(11, 35, 280, 1209)', 1337302)
            LocationData('<areaName>', 'ItemKey(11, 35, 280, 537)', 1337303)
            LocationData('<areaName>', 'ItemKey(11, 35, 280, 905)', 1337304)
            LocationData('<areaName>', 'ItemKey(11, 36, 126, 176)', 1337305)
            LocationData('<areaName>', 'ItemKey(11, 36, 190, 176)', 1337306)
            LocationData('<areaName>', 'ItemKey(11, 36, 254, 176)', 1337307)
            LocationData('<areaName>', 'ItemKey(11, 37, 136, 169)', 1337308)
            LocationData('<areaName>', 'ItemKey(11, 37, 264, 169)', 1337309)
            LocationData('<areaName>', 'ItemKey(11, 4, 136, 169)', 1337310)
            LocationData('<areaName>', 'ItemKey(11, 4, 264, 169)', 1337311)
            LocationData('<areaName>', 'ItemKey(11, 6, 136, 169)', 1337312)
            LocationData('<areaName>', 'ItemKey(11, 6, 264, 169)', 1337313)
            LocationData('<areaName>', 'ItemKey(12, 1, 122, 572)', 1337314)
            LocationData('<areaName>', 'ItemKey(12, 1, 218, 172)', 1337315)
            LocationData('<areaName>', 'ItemKey(12, 1, 234, 572)', 1337316)
            LocationData('<areaName>', 'ItemKey(12, 12, 1034, 140)', 1337317)
            LocationData('<areaName>', 'ItemKey(12, 12, 1466, 140)', 1337318)
            LocationData('<areaName>', 'ItemKey(12, 12, 170, 140)', 1337319)
            LocationData('<areaName>', 'ItemKey(12, 12, 602, 140)', 1337320)
            LocationData('<areaName>', 'ItemKey(12, 16, 122, 140)', 1337321)
            LocationData('<areaName>', 'ItemKey(12, 16, 266, 140)', 1337322)
            LocationData('<areaName>', 'ItemKey(12, 17, 106, 172)', 1337323)
            LocationData('<areaName>', 'ItemKey(12, 17, 282, 172)', 1337324)
            LocationData('<areaName>', 'ItemKey(12, 19, 106, 156)', 1337325)
            LocationData('<areaName>', 'ItemKey(12, 19, 282, 156)', 1337326)
            LocationData('<areaName>', 'ItemKey(12, 21, 122, 892)', 1337327)
            LocationData('<areaName>', 'ItemKey(12, 21, 298, 892)', 1337328)
            LocationData('<areaName>', 'ItemKey(12, 22, 106, 156)', 1337329)
            LocationData('<areaName>', 'ItemKey(12, 22, 282, 156)', 1337330)
            LocationData('<areaName>', 'ItemKey(12, 24, 122, 140)', 1337331)
            LocationData('<areaName>', 'ItemKey(12, 24, 266, 140)', 1337332)
            LocationData('<areaName>', 'ItemKey(12, 25, 122, 140)', 1337333)
            LocationData('<areaName>', 'ItemKey(12, 26, 106, 156)', 1337334)
            LocationData('<areaName>', 'ItemKey(12, 26, 282, 156)', 1337335)
            LocationData('<areaName>', 'ItemKey(12, 3, 138, 124)', 1337336)
            LocationData('<areaName>', 'ItemKey(12, 3, 282, 124)', 1337337)
            LocationData('<areaName>', 'ItemKey(12, 4, 106, 172)', 1337338)
            LocationData('<areaName>', 'ItemKey(12, 4, 106, 572)', 1337339)
            LocationData('<areaName>', 'ItemKey(12, 4, 282, 172)', 1337340)
            LocationData('<areaName>', 'ItemKey(12, 4, 282, 572)', 1337341)
            LocationData('<areaName>', 'ItemKey(12, 5, 122, 156)', 1337342)
            LocationData('<areaName>', 'ItemKey(12, 5, 298, 156)', 1337343)
            LocationData('<areaName>', 'ItemKey(12, 6, 1178, 140)', 1337344)
            LocationData('<areaName>', 'ItemKey(12, 6, 1466, 140)', 1337345)
            LocationData('<areaName>', 'ItemKey(12, 6, 314, 140)', 1337346)
            LocationData('<areaName>', 'ItemKey(12, 6, 602, 140)', 1337347)
            LocationData('<areaName>', 'ItemKey(12, 6, 890, 140)', 1337348)
            LocationData('<areaName>', 'ItemKey(12, 7, 202, 1244)', 1337349)
            LocationData('<areaName>', 'ItemKey(12, 7, 202, 332)', 1337350)
            LocationData('<areaName>', 'ItemKey(12, 9, 122, 892)', 1337351)
            LocationData('<areaName>', 'ItemKey(12, 9, 298, 892)', 1337352)
            LocationData('<areaName>', 'ItemKey(14, 12, 280, 128)', 1337353)
            LocationData('<areaName>', 'ItemKey(14, 12, 520, 128)', 1337354)
            LocationData('<areaName>', 'ItemKey(14, 13, 280, 128)', 1337355)
            LocationData('<areaName>', 'ItemKey(14, 13, 520, 128)', 1337356)
            LocationData('<areaName>', 'ItemKey(14, 15, 280, 128)', 1337357)
            LocationData('<areaName>', 'ItemKey(14, 15, 520, 128)', 1337358)
            LocationData('<areaName>', 'ItemKey(14, 16, 280, 128)', 1337359)
            LocationData('<areaName>', 'ItemKey(14, 16, 520, 128)', 1337360)
            LocationData('<areaName>', 'ItemKey(14, 18, 280, 128)', 1337361)
            LocationData('<areaName>', 'ItemKey(14, 18, 520, 128)', 1337362)
            LocationData('<areaName>', 'ItemKey(14, 19, 280, 128)', 1337363)
            LocationData('<areaName>', 'ItemKey(14, 19, 520, 128)', 1337364)
            LocationData('<areaName>', 'ItemKey(14, 21, 280, 128)', 1337365)
            LocationData('<areaName>', 'ItemKey(14, 21, 520, 128)', 1337366)
            LocationData('<areaName>', 'ItemKey(14, 22, 280, 128)', 1337367)
            LocationData('<areaName>', 'ItemKey(14, 22, 520, 128)', 1337368)
            LocationData('<areaName>', 'ItemKey(15, 0, 1144, 169)', 1337369)
            LocationData('<areaName>', 'ItemKey(15, 0, 1448, 169)', 1337370)
            LocationData('<areaName>', 'ItemKey(15, 0, 344, 169)', 1337371)
            LocationData('<areaName>', 'ItemKey(15, 0, 616, 169)', 1337372)
            LocationData('<areaName>', 'ItemKey(15, 0, 824, 169)', 1337373)
            LocationData('<areaName>', 'ItemKey(15, 0, 88, 169)', 1337374)
            LocationData('<areaName>', 'ItemKey(15, 2, 104, 537)', 1337375)
            LocationData('<areaName>', 'ItemKey(15, 2, 344, 537)', 1337376)
            LocationData('<areaName>', 'ItemKey(16, 0, 152, 128)', 1337377)
            LocationData('<areaName>', 'ItemKey(16, 0, 247, 128)', 1337378)
            LocationData('<areaName>', 'ItemKey(16, 1, 202, 1526)', 1337379)
            LocationData('<areaName>', 'ItemKey(16, 1, 2202, 166)', 1337380)
            LocationData('<areaName>', 'ItemKey(16, 10, 202, 1526)', 1337381)
            LocationData('<areaName>', 'ItemKey(16, 10, 2202, 166)', 1337382)
            LocationData('<areaName>', 'ItemKey(16, 10, 634, 1494)', 1337383)
            LocationData('<areaName>', 'ItemKey(16, 11, 202, 166)', 1337384)
            LocationData('<areaName>', 'ItemKey(16, 11, 2202, 1526)', 1337385)
            LocationData('<areaName>', 'ItemKey(16, 14, 202, 166)', 1337386)
            LocationData('<areaName>', 'ItemKey(16, 15, 154, 886)', 1337387)
            LocationData('<areaName>', 'ItemKey(16, 15, 250, 886)', 1337388)
            LocationData('<areaName>', 'ItemKey(16, 16, 1002, 166)', 1337389)
            LocationData('<areaName>', 'ItemKey(16, 16, 1402, 166)', 1337390)
            LocationData('<areaName>', 'ItemKey(16, 16, 202, 166)', 1337391)
            LocationData('<areaName>', 'ItemKey(16, 16, 602, 166)', 1337392)
            LocationData('<areaName>', 'ItemKey(16, 17, 154, 566)', 1337393)
            LocationData('<areaName>', 'ItemKey(16, 17, 250, 566)', 1337394)
            LocationData('<areaName>', 'ItemKey(16, 18, 1002, 166)', 1337395)
            LocationData('<areaName>', 'ItemKey(16, 18, 1402, 166)', 1337396)
            LocationData('<areaName>', 'ItemKey(16, 18, 202, 166)', 1337397)
            LocationData('<areaName>', 'ItemKey(16, 18, 602, 166)', 1337398)
            LocationData('<areaName>', 'ItemKey(16, 2, 1370, 822)', 1337399)
            LocationData('<areaName>', 'ItemKey(16, 2, 202, 166)', 1337400)
            LocationData('<areaName>', 'ItemKey(16, 2, 2202, 1526)', 1337401)
            LocationData('<areaName>', 'ItemKey(16, 20, 154, 886)', 1337402)
            LocationData('<areaName>', 'ItemKey(16, 20, 250, 886)', 1337403)
            LocationData('<areaName>', 'ItemKey(16, 22, 138, 166)', 1337404)
            LocationData('<areaName>', 'ItemKey(16, 22, 266, 166)', 1337405)
            LocationData('<areaName>', 'ItemKey(16, 23, 154, 566)', 1337406)
            LocationData('<areaName>', 'ItemKey(16, 23, 250, 566)', 1337407)
            LocationData('<areaName>', 'ItemKey(16, 24, 154, 566)', 1337408)
            LocationData('<areaName>', 'ItemKey(16, 24, 250, 566)', 1337409)
            LocationData('<areaName>', 'ItemKey(16, 25, 202, 166)', 1337410)
            LocationData('<areaName>', 'ItemKey(16, 3, 152, 128)', 1337411)
            LocationData('<areaName>', 'ItemKey(16, 3, 247, 128)', 1337412)
            LocationData('<areaName>', 'ItemKey(16, 5, 152, 128)', 1337413)
            LocationData('<areaName>', 'ItemKey(16, 5, 247, 128)', 1337414)
            LocationData('<areaName>', 'ItemKey(16, 6, 154, 2486)', 1337415)
            LocationData('<areaName>', 'ItemKey(16, 6, 250, 2486)', 1337416)
            LocationData('<areaName>', 'ItemKey(16, 7, 202, 166)', 1337417)
            LocationData('<areaName>', 'ItemKey(16, 7, 602, 166)', 1337418)
            LocationData('<areaName>', 'ItemKey(16, 8, 202, 166)', 1337419)
            LocationData('<areaName>', 'ItemKey(16, 9, 202, 166)', 1337420)
            LocationData('<areaName>', 'ItemKey(16, 9, 602, 166)', 1337421)
            LocationData('<areaName>', 'ItemKey(2, 0, 1432, 116)', 1337422)
            LocationData('<areaName>', 'ItemKey(2, 0, 504, 116)', 1337423)
            LocationData('<areaName>', 'ItemKey(2, 1, 1208, 116)', 1337424)
            LocationData('<areaName>', 'ItemKey(2, 1, 1320, 116)', 1337425)
            LocationData('<areaName>', 'ItemKey(2, 1, 280, 116)', 1337426)
            LocationData('<areaName>', 'ItemKey(2, 1, 392, 116)', 1337427)
            LocationData('<areaName>', 'ItemKey(2, 1, 616, 116)', 1337428)
            LocationData('<areaName>', 'ItemKey(2, 1, 984, 116)', 1337429)
            LocationData('<areaName>', 'ItemKey(2, 10, 2536, 116)', 1337430)
            LocationData('<areaName>', 'ItemKey(2, 10, 664, 116)', 1337431)
            LocationData('<areaName>', 'ItemKey(2, 12, 168, 256)', 1337432)
            LocationData('<areaName>', 'ItemKey(2, 12, 632, 256)', 1337433)
            LocationData('<areaName>', 'ItemKey(2, 16, 200, 192)', 1337434)
            LocationData('<areaName>', 'ItemKey(2, 17, 168, 192)', 1337435)
            LocationData('<areaName>', 'ItemKey(2, 17, 648, 192)', 1337436)
            LocationData('<areaName>', 'ItemKey(2, 19, 1048, 192)', 1337437)
            LocationData('<areaName>', 'ItemKey(2, 19, 168, 192)', 1337438)
            LocationData('<areaName>', 'ItemKey(2, 20, 168, 1200)', 1337439)
            LocationData('<areaName>', 'ItemKey(2, 20, 168, 1840)', 1337440)
            LocationData('<areaName>', 'ItemKey(2, 20, 168, 208)', 1337441)
            LocationData('<areaName>', 'ItemKey(2, 20, 168, 2160)', 1337442)
            LocationData('<areaName>', 'ItemKey(2, 20, 168, 2512)', 1337443)
            LocationData('<areaName>', 'ItemKey(2, 20, 168, 880)', 1337444)
            LocationData('<areaName>', 'ItemKey(2, 23, 200, 560)', 1337445)
            LocationData('<areaName>', 'ItemKey(2, 23, 648, 560)', 1337446)
            LocationData('<areaName>', 'ItemKey(2, 23, 984, 560)', 1337447)
            LocationData('<areaName>', 'ItemKey(2, 24, 1048, 192)', 1337448)
            LocationData('<areaName>', 'ItemKey(2, 24, 88, 192)', 1337449)
            LocationData('<areaName>', 'ItemKey(2, 31, 208, 96)', 1337450)
            LocationData('<areaName>', 'ItemKey(2, 32, 120, 192)', 1337451)
            LocationData('<areaName>', 'ItemKey(2, 34, 232, 1520)', 1337452)
            LocationData('<areaName>', 'ItemKey(2, 34, 232, 2512)', 1337453)
            LocationData('<areaName>', 'ItemKey(2, 36, 312, 208)', 1337454)
            LocationData('<areaName>', 'ItemKey(2, 37, 152, 192)', 1337455)
            LocationData('<areaName>', 'ItemKey(2, 37, 728, 192)', 1337456)
            LocationData('<areaName>', 'ItemKey(2, 39, 264, 192)', 1337457)
            LocationData('<areaName>', 'ItemKey(2, 39, 552, 192)', 1337458)
            LocationData('<areaName>', 'ItemKey(2, 4, 112, 512)', 1337459)
            LocationData('<areaName>', 'ItemKey(2, 4, 112, 96)', 1337460)
            LocationData('<areaName>', 'ItemKey(2, 4, 288, 512)', 1337461)
            LocationData('<areaName>', 'ItemKey(2, 4, 288, 96)', 1337462)
            LocationData('<areaName>', 'ItemKey(2, 41, 176, 208)', 1337463)
            LocationData('<areaName>', 'ItemKey(2, 41, 624, 208)', 1337464)
            LocationData('<areaName>', 'ItemKey(2, 42, 112, 512)', 1337465)
            LocationData('<areaName>', 'ItemKey(2, 42, 112, 96)', 1337466)
            LocationData('<areaName>', 'ItemKey(2, 42, 288, 512)', 1337467)
            LocationData('<areaName>', 'ItemKey(2, 42, 288, 96)', 1337468)
            LocationData('<areaName>', 'ItemKey(2, 43, 200, 192)', 1337469)
            LocationData('<areaName>', 'ItemKey(2, 43, 600, 192)', 1337470)
            LocationData('<areaName>', 'ItemKey(2, 44, 1072, 496)', 1337471)
            LocationData('<areaName>', 'ItemKey(2, 44, 112, 496)', 1337472)
            LocationData('<areaName>', 'ItemKey(2, 44, 400, 496)', 1337473)
            LocationData('<areaName>', 'ItemKey(2, 44, 400, 96)', 1337474)
            LocationData('<areaName>', 'ItemKey(2, 44, 624, 496)', 1337475)
            LocationData('<areaName>', 'ItemKey(2, 44, 736, 96)', 1337476)
            LocationData('<areaName>', 'ItemKey(2, 44, 848, 496)', 1337477)
            LocationData('<areaName>', 'ItemKey(2, 47, 128, 96)', 1337478)
            LocationData('<areaName>', 'ItemKey(2, 47, 240, 96)', 1337479)
            LocationData('<areaName>', 'ItemKey(2, 5, 120, 100)', 1337480)
            LocationData('<areaName>', 'ItemKey(2, 5, 280, 100)', 1337481)
            LocationData('<areaName>', 'ItemKey(2, 51, 128, 96)', 1337482)
            LocationData('<areaName>', 'ItemKey(2, 51, 240, 96)', 1337483)
            LocationData('<areaName>', 'ItemKey(2, 56, 1072, 96)', 1337484)
            LocationData('<areaName>', 'ItemKey(2, 56, 176, 96)', 1337485)
            LocationData('<areaName>', 'ItemKey(2, 56, 400, 96)', 1337486)
            LocationData('<areaName>', 'ItemKey(2, 56, 624, 96)', 1337487)
            LocationData('<areaName>', 'ItemKey(2, 56, 848, 96)', 1337488)
            LocationData('<areaName>', 'ItemKey(2, 57, 200, 100)', 1337489)
            LocationData('<areaName>', 'ItemKey(2, 58, 160, 96)', 1337490)
            LocationData('<areaName>', 'ItemKey(2, 58, 272, 96)', 1337491)
            LocationData('<areaName>', 'ItemKey(2, 59, 112, 96)', 1337492)
            LocationData('<areaName>', 'ItemKey(2, 59, 288, 96)', 1337493)
            LocationData('<areaName>', 'ItemKey(2, 60, 184, 100)', 1337494)
            LocationData('<areaName>', 'ItemKey(2, 9, 200, 192)', 1337495)
            LocationData('<areaName>', 'ItemKey(2, 9, 600, 192)', 1337496)
            LocationData('<areaName>', 'ItemKey(3, 1, 235, 125)', 1337497)
            LocationData('<areaName>', 'ItemKey(3, 12, 1019, 109)', 1337498)
            LocationData('<areaName>', 'ItemKey(3, 13, 1227, 125)', 1337499)
            LocationData('<areaName>', 'ItemKey(3, 13, 251, 125)', 1337500)
            LocationData('<areaName>', 'ItemKey(3, 13, 747, 125)', 1337501)
            LocationData('<areaName>', 'ItemKey(3, 2, 235, 157)', 1337502)
            LocationData('<areaName>', 'ItemKey(3, 4, 1691, 77)', 1337503)
            LocationData('<areaName>', 'ItemKey(3, 4, 411, 77)', 1337504)
            LocationData('<areaName>', 'ItemKey(3, 7, 1307, 125)', 1337505)
            LocationData('<areaName>', 'ItemKey(11, 17, 120, 505)', 1337506)
            LocationData('<areaName>', 'ItemKey(4, 10, 122, 161)', 1337507)
            LocationData('<areaName>', 'ItemKey(4, 10, 266, 161)', 1337508)
            LocationData('<areaName>', 'ItemKey(4, 11, 138, 161)', 1337509)
            LocationData('<areaName>', 'ItemKey(4, 11, 282, 161)', 1337510)
            LocationData('<areaName>', 'ItemKey(6, 23, 519, 73)', 1337511)
            LocationData('<areaName>', 'ItemKey(6, 23, 599, 73)', 1337512)
            LocationData('<areaName>', 'ItemKey(4, 22, 298, 161)', 1337513)
            LocationData('<areaName>', 'ItemKey(4, 3, 234, 561)', 1337514)
            LocationData('<areaName>', 'ItemKey(4, 3, 570, 561)', 1337515)
            LocationData('<areaName>', 'ItemKey(4, 4, 1418, 129)', 1337516)
            LocationData('<areaName>', 'ItemKey(4, 4, 378, 129)', 1337517)
            LocationData('<areaName>', 'ItemKey(4, 4, 922, 129)', 1337518)
            LocationData('<areaName>', 'ItemKey(4, 5, 202, 497)', 1337519)
            LocationData('<areaName>', 'ItemKey(4, 5, 312, 264)', 1337520)
            LocationData('<areaName>', 'ItemKey(4, 5, 88, 264)', 1337521)
            LocationData('<areaName>', 'ItemKey(4, 6, 234, 561)', 1337522)
            LocationData('<areaName>', 'ItemKey(4, 6, 280, 104)', 1337523)
            LocationData('<areaName>', 'ItemKey(4, 6, 520, 104)', 1337524)
            LocationData('<areaName>', 'ItemKey(4, 6, 570, 561)', 1337525)
            LocationData('<areaName>', 'ItemKey(4, 7, 202, 497)', 1337526)
            LocationData('<areaName>', 'ItemKey(4, 7, 312, 264)', 1337527)
            LocationData('<areaName>', 'ItemKey(4, 7, 88, 264)', 1337528)
            LocationData('<areaName>', 'ItemKey(4, 8, 1418, 129)', 1337529)
            LocationData('<areaName>', 'ItemKey(4, 8, 362, 129)', 1337530)
            LocationData('<areaName>', 'ItemKey(4, 8, 858, 129)', 1337531)
            LocationData('<areaName>', 'ItemKey(5, 0, 1032, 173)', 1337532)
            LocationData('<areaName>', 'ItemKey(5, 0, 184, 173)', 1337533)
            LocationData('<areaName>', 'ItemKey(5, 0, 456, 173)', 1337534)
            LocationData('<areaName>', 'ItemKey(5, 0, 744, 173)', 1337535)
            LocationData('<areaName>', 'ItemKey(5, 1, 200, 188)', 1337536)
            LocationData('<areaName>', 'ItemKey(5, 1, 200, 444)', 1337537)
            LocationData('<areaName>', 'ItemKey(5, 1, 200, 700)', 1337538)
            LocationData('<areaName>', 'ItemKey(5, 11, 122, 161)', 1337539)
            LocationData('<areaName>', 'ItemKey(5, 11, 298, 161)', 1337540)
            LocationData('<areaName>', 'ItemKey(5, 13, 200, 188)', 1337541)
            LocationData('<areaName>', 'ItemKey(5, 13, 200, 444)', 1337542)
            LocationData('<areaName>', 'ItemKey(5, 16, 298, 129)', 1337543)
            LocationData('<areaName>', 'ItemKey(5, 16, 490, 129)', 1337544)
            LocationData('<areaName>', 'ItemKey(5, 18, 200, 172)', 1337545)
            LocationData('<areaName>', 'ItemKey(5, 18, 200, 492)', 1337546)
            LocationData('<areaName>', 'ItemKey(5, 18, 200, 748)', 1337547)
            LocationData('<areaName>', 'ItemKey(5, 19, 168, 157)', 1337548)
            LocationData('<areaName>', 'ItemKey(5, 19, 424, 157)', 1337549)
            LocationData('<areaName>', 'ItemKey(5, 19, 728, 157)', 1337550)
            LocationData('<areaName>', 'ItemKey(5, 19, 984, 157)', 1337551)
            LocationData('<areaName>', 'ItemKey(5, 21, 184, 157)', 1337552)
            LocationData('<areaName>', 'ItemKey(5, 21, 440, 157)', 1337553)
            LocationData('<areaName>', 'ItemKey(5, 21, 728, 157)', 1337554)
            LocationData('<areaName>', 'ItemKey(5, 21, 984, 157)', 1337555)
            LocationData('<areaName>', 'ItemKey(5, 22, 152, 141)', 1337556)
            LocationData('<areaName>', 'ItemKey(5, 22, 264, 141)', 1337557)
            LocationData('<areaName>', 'ItemKey(5, 29, 200, 204)', 1337558)
            LocationData('<areaName>', 'ItemKey(5, 29, 200, 460)', 1337559)
            LocationData('<areaName>', 'ItemKey(5, 3, 1000, 445)', 1337560)
            LocationData('<areaName>', 'ItemKey(5, 3, 1016, 205)', 1337561)
            LocationData('<areaName>', 'ItemKey(5, 3, 184, 205)', 1337562)
            LocationData('<areaName>', 'ItemKey(5, 3, 200, 445)', 1337563)
            LocationData('<areaName>', 'ItemKey(5, 3, 552, 285)', 1337564)
            LocationData('<areaName>', 'ItemKey(5, 3, 648, 285)', 1337565)
            LocationData('<areaName>', 'ItemKey(5, 32, 488, 141)', 1337566)
            LocationData('<areaName>', 'ItemKey(5, 33, 122, 161)', 1337567)
            LocationData('<areaName>', 'ItemKey(5, 33, 298, 161)', 1337568)
            LocationData('<areaName>', 'ItemKey(5, 35, 298, 129)', 1337569)
            LocationData('<areaName>', 'ItemKey(5, 35, 490, 129)', 1337570)
            LocationData('<areaName>', 'ItemKey(5, 39, 200, 157)', 1337571)
            LocationData('<areaName>', 'ItemKey(5, 39, 312, 157)', 1337572)
            LocationData('<areaName>', 'ItemKey(5, 39, 488, 157)', 1337573)
            LocationData('<areaName>', 'ItemKey(5, 39, 600, 157)', 1337574)
            LocationData('<areaName>', 'ItemKey(5, 4, 168, 173)', 1337575)
            LocationData('<areaName>', 'ItemKey(5, 4, 232, 173)', 1337576)
            LocationData('<areaName>', 'ItemKey(5, 41, 152, 173)', 1337577)
            LocationData('<areaName>', 'ItemKey(5, 41, 232, 173)', 1337578)
            LocationData('<areaName>', 'ItemKey(5, 43, 200, 172)', 1337579)
            LocationData('<areaName>', 'ItemKey(5, 43, 200, 428)', 1337580)
            LocationData('<areaName>', 'ItemKey(5, 44, 122, 161)', 1337581)
            LocationData('<areaName>', 'ItemKey(5, 44, 298, 161)', 1337582)
            LocationData('<areaName>', 'ItemKey(5, 6, 264, 141)', 1337583)
            LocationData('<areaName>', 'ItemKey(5, 6, 552, 141)', 1337584)
            LocationData('<areaName>', 'ItemKey(5, 6, 808, 141)', 1337585)
            LocationData('<areaName>', 'ItemKey(5, 7, 200, 188)', 1337586)
            LocationData('<areaName>', 'ItemKey(5, 7, 200, 444)', 1337587)
            LocationData('<areaName>', 'ItemKey(6, 1, 1114, 257)', 1337588)
            LocationData('<areaName>', 'ItemKey(6, 1, 1402, 257)', 1337589)
            LocationData('<areaName>', 'ItemKey(6, 1, 250, 257)', 1337590)
            LocationData('<areaName>', 'ItemKey(6, 1, 538, 257)', 1337591)
            LocationData('<areaName>', 'ItemKey(6, 1, 826, 257)', 1337592)
            LocationData('<areaName>', 'ItemKey(6, 10, 138, 253)', 1337593)
            LocationData('<areaName>', 'ItemKey(6, 10, 266, 253)', 1337594)
            LocationData('<areaName>', 'ItemKey(6, 11, 298, 173)', 1337595)
            LocationData('<areaName>', 'ItemKey(6, 16, 202, 157)', 1337596)
            LocationData('<areaName>', 'ItemKey(6, 16, 346, 157)', 1337597)
            LocationData('<areaName>', 'ItemKey(6, 16, 58, 157)', 1337598)
            LocationData('<areaName>', 'ItemKey(6, 2, 74, 349)', 1337599)
            LocationData('<areaName>', 'ItemKey(6, 2, 74, 541)', 1337600)
            LocationData('<areaName>', 'ItemKey(6, 21, 154, 893)', 1337601)
            LocationData('<areaName>', 'ItemKey(6, 21, 250, 893)', 1337602)
            LocationData('<areaName>', 'ItemKey(6, 22, 154, 173)', 1337603)
            LocationData('<areaName>', 'ItemKey(6, 22, 250, 173)', 1337604)
            LocationData('<areaName>', 'ItemKey(6, 24, 122, 141)', 1337605)
            LocationData('<areaName>', 'ItemKey(6, 24, 266, 141)', 1337606)
            LocationData('<areaName>', 'ItemKey(6, 25, 135, 73)', 1337607)
            LocationData('<areaName>', 'ItemKey(6, 25, 87, 73)', 1337608)
            LocationData('<areaName>', 'ItemKey(6, 27, 394, 157)', 1337609)
            LocationData('<areaName>', 'ItemKey(6, 27, 394, 541)', 1337610)
            LocationData('<areaName>', 'ItemKey(6, 3, 330, 173)', 1337611)
            LocationData('<areaName>', 'ItemKey(6, 7, 202, 1245)', 1337612)
            LocationData('<areaName>', 'ItemKey(6, 7, 202, 333)', 1337613)
            LocationData('<areaName>', 'ItemKey(6, 9, 199, 217)', 1337614)
            LocationData('<areaName>', 'ItemKey(6, 9, 199, 537)', 1337615)
            LocationData('<areaName>', 'ItemKey(6, 9, 199, 841)', 1337616)
            LocationData('<areaName>', 'ItemKey(7, 1, 2568, 107)', 1337617)
            LocationData('<areaName>', 'ItemKey(7, 11, 184, 160)', 1337618)
            LocationData('<areaName>', 'ItemKey(7, 11, 248, 160)', 1337619)
            LocationData('<areaName>', 'ItemKey(7, 11, 312, 160)', 1337620)
            LocationData('<areaName>', 'ItemKey(7, 14, 104, 283)', 1337621)
            LocationData('<areaName>', 'ItemKey(7, 15, 200, 187)', 1337622)
            LocationData('<areaName>', 'ItemKey(7, 15, 2104, 187)', 1337623)
            LocationData('<areaName>', 'ItemKey(7, 15, 2200, 187)', 1337624)
            LocationData('<areaName>', 'ItemKey(7, 15, 296, 187)', 1337625)
            LocationData('<areaName>', 'ItemKey(7, 17, 312, 267)', 1337626)
            LocationData('<areaName>', 'ItemKey(7, 2, 1480, 688)', 1337627)
            LocationData('<areaName>', 'ItemKey(7, 2, 152, 688)', 1337628)
            LocationData('<areaName>', 'ItemKey(7, 2, 936, 656)', 1337629)
            LocationData('<areaName>', 'ItemKey(7, 20, 184, 240)', 1337630)
            LocationData('<areaName>', 'ItemKey(7, 3, 120, 379)', 1337631)
            LocationData('<areaName>', 'ItemKey(7, 3, 328, 736)', 1337632)
            LocationData('<areaName>', 'ItemKey(7, 3, 392, 1104)', 1337633)
            LocationData('<areaName>', 'ItemKey(7, 4, 328, 155)', 1337634)
            LocationData('<areaName>', 'ItemKey(7, 6, 1336, 256)', 1337635)
            LocationData('<areaName>', 'ItemKey(7, 6, 472, 432)', 1337636)
            LocationData('<areaName>', 'ItemKey(7, 6, 568, 192)', 1337637)
            LocationData('<areaName>', 'ItemKey(7, 7, 1112, 592)', 1337638)
            LocationData('<areaName>', 'ItemKey(7, 7, 1544, 592)', 1337639)
            LocationData('<areaName>', 'ItemKey(7, 7, 1704, 155)', 1337640)
            LocationData('<areaName>', 'ItemKey(7, 7, 488, 576)', 1337641)
            LocationData('<areaName>', 'ItemKey(7, 7, 824, 608)', 1337642)
            LocationData('<areaName>', 'ItemKey(8, 0, 1032, 140)', 1337643)
            LocationData('<areaName>', 'ItemKey(8, 0, 1176, 140)', 1337644)
            LocationData('<areaName>', 'ItemKey(8, 0, 1416, 140)', 1337645)
            LocationData('<areaName>', 'ItemKey(8, 0, 1576, 140)', 1337646)
            LocationData('<areaName>', 'ItemKey(8, 0, 424, 140)', 1337647)
            LocationData('<areaName>', 'ItemKey(8, 0, 568, 140)', 1337648)
            LocationData('<areaName>', 'ItemKey(7, 7, 184, 123)', 1337649)
            LocationData('<areaName>', 'ItemKey(8, 10, 120, 124)', 1337650)
            LocationData('<areaName>', 'ItemKey(8, 10, 184, 828)', 1337651)
            LocationData('<areaName>', 'ItemKey(8, 10, 232, 828)', 1337652)
            LocationData('<areaName>', 'ItemKey(8, 10, 296, 124)', 1337653)
            LocationData('<areaName>', 'ItemKey(8, 11, 153, 140)', 1337654)
            LocationData('<areaName>', 'ItemKey(8, 11, 409, 140)', 1337655)
            LocationData('<areaName>', 'ItemKey(8, 11, 793, 140)', 1337656)
            LocationData('<areaName>', 'ItemKey(8, 11, 985, 140)', 1337657)
            LocationData('<areaName>', 'ItemKey(8, 14, 1050, 157)', 1337658)
            LocationData('<areaName>', 'ItemKey(8, 14, 1258, 157)', 1337659)
            LocationData('<areaName>', 'ItemKey(8, 14, 1466, 157)', 1337660)
            LocationData('<areaName>', 'ItemKey(8, 14, 218, 157)', 1337661)
            LocationData('<areaName>', 'ItemKey(8, 14, 426, 157)', 1337662)
            LocationData('<areaName>', 'ItemKey(8, 14, 634, 157)', 1337663)
            LocationData('<areaName>', 'ItemKey(8, 16, 104, 140)', 1337664)
            LocationData('<areaName>', 'ItemKey(8, 16, 104, 1756)', 1337665)
            LocationData('<areaName>', 'ItemKey(8, 16, 1096, 140)', 1337666)
            LocationData('<areaName>', 'ItemKey(8, 16, 1096, 1756)', 1337667)
            LocationData('<areaName>', 'ItemKey(8, 16, 1096, 780)', 1337668)
            LocationData('<areaName>', 'ItemKey(8, 16, 920, 780)', 1337669)
            LocationData('<areaName>', 'ItemKey(8, 19, 520, 876)', 1337670)
            LocationData('<areaName>', 'ItemKey(8, 19, 536, 140)', 1337671)
            LocationData('<areaName>', 'ItemKey(8, 19, 696, 876)', 1337672)
            LocationData('<areaName>', 'ItemKey(8, 19, 712, 140)', 1337673)
            LocationData('<areaName>', 'ItemKey(8, 20, 1001, 140)', 1337674)
            LocationData('<areaName>', 'ItemKey(8, 20, 297, 140)', 1337675)
            LocationData('<areaName>', 'ItemKey(8, 20, 617, 140)', 1337676)
            LocationData('<areaName>', 'ItemKey(8, 31, 105, 540)', 1337677)
            LocationData('<areaName>', 'ItemKey(8, 31, 265, 572)', 1337678)
            LocationData('<areaName>', 'ItemKey(8, 31, 281, 1516)', 1337679)
            LocationData('<areaName>', 'ItemKey(8, 31, 329, 124)', 1337680)
            LocationData('<areaName>', 'ItemKey(8, 31, 73, 124)', 1337681)
            LocationData('<areaName>', 'ItemKey(8, 37, 1032, 124)', 1337682)
            LocationData('<areaName>', 'ItemKey(8, 37, 1176, 124)', 1337683)
            LocationData('<areaName>', 'ItemKey(8, 37, 424, 124)', 1337684)
            LocationData('<areaName>', 'ItemKey(8, 37, 584, 124)', 1337685)
            LocationData('<areaName>', 'ItemKey(8, 38, 120, 124)', 1337686)
            LocationData('<areaName>', 'ItemKey(8, 38, 120, 508)', 1337687)
            LocationData('<areaName>', 'ItemKey(8, 38, 168, 828)', 1337688)
            LocationData('<areaName>', 'ItemKey(8, 38, 248, 828)', 1337689)
            LocationData('<areaName>', 'ItemKey(8, 38, 296, 124)', 1337690)
            LocationData('<areaName>', 'ItemKey(8, 43, 1050, 157)', 1337691)
            LocationData('<areaName>', 'ItemKey(8, 43, 1258, 157)', 1337692)
            LocationData('<areaName>', 'ItemKey(8, 43, 1466, 157)', 1337693)
            LocationData('<areaName>', 'ItemKey(8, 43, 218, 157)', 1337694)
            LocationData('<areaName>', 'ItemKey(8, 43, 426, 157)', 1337695)
            LocationData('<areaName>', 'ItemKey(8, 46, 104, 140)', 1337696)
            LocationData('<areaName>', 'ItemKey(8, 46, 296, 140)', 1337697)
            LocationData('<areaName>', 'ItemKey(8, 47, 122, 157)', 1337698)
            LocationData('<areaName>', 'ItemKey(8, 47, 266, 157)', 1337699)
            LocationData('<areaName>', 'ItemKey(8, 5, 121, 156)', 1337700)
            LocationData('<areaName>', 'ItemKey(8, 5, 1241, 156)', 1337701)
            LocationData('<areaName>', 'ItemKey(8, 5, 377, 156)', 1337702)
            LocationData('<areaName>', 'ItemKey(8, 5, 985, 156)', 1337703)
            LocationData('<areaName>', 'ItemKey(8, 8, 282, 157)', 1337704)
            LocationData('<areaName>', 'ItemKey(8, 8, 392, 348)', 1337705)
            LocationData('<areaName>', 'ItemKey(8, 8, 392, 540)', 1337706)
            LocationData('<areaName>', 'ItemKey(8, 8, 506, 157)', 1337707)
            LocationData('<areaName>', 'ItemKey(8, 9, 104, 124)', 1337708)
            LocationData('<areaName>', 'ItemKey(8, 9, 104, 508)', 1337709)
            LocationData('<areaName>', 'ItemKey(8, 9, 280, 124)', 1337710)
            LocationData('<areaName>', 'ItemKey(8, 9, 280, 508)', 1337711)
            LocationData('<areaName>', 'ItemKey(9, 0, 1032, 140)', 1337712)
            LocationData('<areaName>', 'ItemKey(9, 0, 1176, 140)', 1337713)
            LocationData('<areaName>', 'ItemKey(9, 0, 1416, 140)', 1337714)
            LocationData('<areaName>', 'ItemKey(9, 0, 1576, 140)', 1337715)
            LocationData('<areaName>', 'ItemKey(9, 0, 424, 140)', 1337716)
            LocationData('<areaName>', 'ItemKey(9, 0, 568, 140)', 1337717)
            LocationData('<areaName>', 'ItemKey(9, 10, 120, 124)', 1337718)
            LocationData('<areaName>', 'ItemKey(9, 10, 184, 828)', 1337719)
            LocationData('<areaName>', 'ItemKey(9, 10, 296, 124)', 1337720)
            LocationData('<areaName>', 'ItemKey(9, 11, 153, 140)', 1337721)
            LocationData('<areaName>', 'ItemKey(9, 11, 409, 140)', 1337722)
            LocationData('<areaName>', 'ItemKey(9, 11, 793, 140)', 1337723)
            LocationData('<areaName>', 'ItemKey(9, 11, 985, 140)', 1337724)
            LocationData('<areaName>', 'ItemKey(9, 14, 1050, 157)', 1337725)
            LocationData('<areaName>', 'ItemKey(9, 14, 1258, 157)', 1337726)
            LocationData('<areaName>', 'ItemKey(9, 14, 1466, 157)', 1337727)
            LocationData('<areaName>', 'ItemKey(9, 14, 218, 157)', 1337728)
            LocationData('<areaName>', 'ItemKey(9, 14, 426, 157)', 1337729)
            LocationData('<areaName>', 'ItemKey(9, 14, 634, 157)', 1337730)
            LocationData('<areaName>', 'ItemKey(9, 16, 104, 140)', 1337731)
            LocationData('<areaName>', 'ItemKey(9, 16, 104, 1756)', 1337732)
            LocationData('<areaName>', 'ItemKey(9, 16, 1096, 140)', 1337733)
            LocationData('<areaName>', 'ItemKey(9, 16, 1096, 1756)', 1337734)
            LocationData('<areaName>', 'ItemKey(9, 16, 1096, 780)', 1337735)
            LocationData('<areaName>', 'ItemKey(9, 16, 920, 780)', 1337736)
            LocationData('<areaName>', 'ItemKey(8, 43, 634, 157)', 1337737)
            LocationData('<areaName>', 'ItemKey(9, 19, 520, 876)', 1337738)
            LocationData('<areaName>', 'ItemKey(9, 19, 536, 140)', 1337739)
            LocationData('<areaName>', 'ItemKey(9, 19, 696, 876)', 1337740)
            LocationData('<areaName>', 'ItemKey(9, 20, 1001, 140)', 1337741)
            LocationData('<areaName>', 'ItemKey(9, 20, 297, 140)', 1337742)
            LocationData('<areaName>', 'ItemKey(9, 20, 617, 140)', 1337743)
            LocationData('<areaName>', 'ItemKey(9, 31, 329, 124)', 1337744)
            LocationData('<areaName>', 'ItemKey(9, 31, 73, 124)', 1337745)
            LocationData('<areaName>', 'ItemKey(9, 37, 1032, 124)', 1337746)
            LocationData('<areaName>', 'ItemKey(9, 37, 1176, 124)', 1337747)
            LocationData('<areaName>', 'ItemKey(9, 37, 424, 124)', 1337748)
            LocationData('<areaName>', 'ItemKey(9, 37, 584, 124)', 1337749)
            LocationData('<areaName>', 'ItemKey(9, 38, 120, 124)', 1337750)
            LocationData('<areaName>', 'ItemKey(9, 38, 120, 508)', 1337751)
            LocationData('<areaName>', 'ItemKey(9, 38, 168, 828)', 1337752)
            LocationData('<areaName>', 'ItemKey(9, 38, 248, 828)', 1337753)
            LocationData('<areaName>', 'ItemKey(9, 38, 296, 124)', 1337754)
            LocationData('<areaName>', 'ItemKey(9, 43, 1050, 157)', 1337755)
            LocationData('<areaName>', 'ItemKey(9, 43, 1258, 157)', 1337756)
            LocationData('<areaName>', 'ItemKey(9, 43, 1466, 157)', 1337757)
            LocationData('<areaName>', 'ItemKey(9, 43, 218, 157)', 1337758)
            LocationData('<areaName>', 'ItemKey(9, 43, 426, 157)', 1337759)
            LocationData('<areaName>', 'ItemKey(9, 43, 634, 157)', 1337760)
            LocationData('<areaName>', 'ItemKey(9, 46, 104, 140)', 1337761)
            LocationData('<areaName>', 'ItemKey(9, 46, 296, 140)', 1337762)
            LocationData('<areaName>', 'ItemKey(9, 47, 122, 157)', 1337763)
            LocationData('<areaName>', 'ItemKey(9, 47, 266, 157)', 1337764)
            LocationData('<areaName>', 'ItemKey(9, 8, 282, 157)', 1337765)
            LocationData('<areaName>', 'ItemKey(9, 8, 392, 348)', 1337766)
            LocationData('<areaName>', 'ItemKey(9, 8, 392, 540)', 1337767)
            LocationData('<areaName>', 'ItemKey(9, 8, 506, 157)', 1337768)
            LocationData('<areaName>', 'ItemKey(9, 9, 104, 124)', 1337769)
            LocationData('<areaName>', 'ItemKey(9, 9, 104, 508)', 1337770)
            LocationData('<areaName>', 'ItemKey(9, 9, 280, 124)', 1337771)
            LocationData('<areaName>', 'ItemKey(9, 9, 280, 508)', 1337772)
            LocationData('<areaName>', 'ItemKey(1, 7, 472, 95)', 1337773)
            LocationData('<areaName>', 'ItemKey(1, 7, 664, 95)', 1337774)
            LocationData('<areaName>', 'ItemKey(1, 7, 856, 95)', 1337775)
            LocationData('<areaName>', 'ItemKey(1, 7, 1048, 95)', 1337776)
            LocationData('<areaName>', 'ItemKey(1, 7, 1240, 95)', 1337777)
            LocationData('<areaName>', 'ItemKey(1, 7, 1432, 95)', 1337778)
            LocationData('<areaName>', 'ItemKey(1, 7, 1624, 95)', 1337779)
            LocationData('<areaName>', 'ItemKey(6, 23, 679, 73)', 1337780)
            LocationData('<areaName>', 'ItemKey(9, 19, 712, 140)', 1337781)
        )

    # Reserved through 1337999
 
    return location_table
