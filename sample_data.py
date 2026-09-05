"""Sample, background and detector metadata used by the project."""

SAMPLE_INFO = {
    "Barrowsoil.TKA": {
        "detector": 1,
        "sample_mass_g": 648.2,
        "expected_live_time_s": 574951,
        "filled_height_cm": 2.6,
    },
    "Barrowswede.TKA": {
        "detector": 1,
        "sample_mass_g": 245.6,
        "expected_live_time_s": 9266,
        "filled_height_cm": 2.4,
    },
    "CARROTSDET1.TKA": {
        "detector": 1,
        "sample_mass_g": 474.0,
        "expected_live_time_s": 575055,
        "filled_height_cm": 2.4,
    },
    "KENDAL CABBAGE DET1.TKA": {
        "detector": 1,
        "sample_mass_g": 362.5,
        "expected_live_time_s": 12880,
        "filled_height_cm": 2.7,
    },
    "WALES LEEK DET1.TKA": {
        "detector": 1,
        "sample_mass_g": 278.0,
        "expected_live_time_s": 12673,
        "filled_height_cm": 2.4,
    },
    "Walessoil.TKA": {
        "detector": 1,
        "sample_mass_g": 508.2,
        "expected_live_time_s": 10481,
        "filled_height_cm": 2.7,
    },
    "WAPLEY ONION DET1.TKA": {
        "detector": 1,
        "sample_mass_g": 453.1,
        "expected_live_time_s": 575097,
        "filled_height_cm": 2.7,
    },
    "WAPLEYSOIL DET1.TKA": {
        "detector": 1,
        "sample_mass_g": 1228.2,
        "expected_live_time_s": 12580,
        "filled_height_cm": 2.7,
    },
    "HPGe2bathcarrots.TKA": {
        "detector": 2,
        "sample_mass_g": 92.8,
        "expected_live_time_s": 10492,
        "filled_height_cm": 1.2,
    },
    "HPGe2bathpotato.TKA": {
        "detector": 2,
        "sample_mass_g": 175.8,
        "expected_live_time_s": 9271,
        "filled_height_cm": 1.2,
    },
    "hpge2bwsoil.TKA": {
        "detector": 2,
        "sample_mass_g": 1358.2,
        "expected_live_time_s": 575449,
        "filled_height_cm": 2.7,
    },
    "hpge2essex.TKA": {
        "detector": 2,
        "sample_mass_g": 1068.2,
        "expected_live_time_s": 13189,
        "filled_height_cm": 2.7,
    },
    "hpge2kendalsoil.TKA": {
        "detector": 2,
        "sample_mass_g": 1158.2,
        "expected_live_time_s": 12888,
        "filled_height_cm": 2.7,
    },
    "hpge2mankypotatoes.TKA": {
        "detector": 2,
        "sample_mass_g": 658.2,
        "expected_live_time_s": 575559,
        "filled_height_cm": 2.7,
    },
    "hpge2walescarrot.TKA": {
        "detector": 2,
        "sample_mass_g": 274.7,
        "expected_live_time_s": 12684,
        "filled_height_cm": 2.7,
    },
    "hpge2essexcarrots.TKA": {
        "detector": 2,
        "sample_mass_g": 401.7,
        "expected_live_time_s": 575562,
        "filled_height_cm": 2.5,
    },
    "BARROWLEEK.TKA": {
        "detector": 3,
        "sample_mass_g": 238.5,
        "expected_live_time_s": 12666,
        "filled_height_cm": 3.9,
    },
    "BATHSOILDET3.TKA": {
        "detector": 3,
        "sample_mass_g": 76.2,
        "expected_live_time_s": 10477,
        "filled_height_cm": 1.7,
    },
    "BathSWTPOT.TKA": {
        "detector": 3,
        "sample_mass_g": 140.2,
        "expected_live_time_s": 9260,
        "filled_height_cm": 2.7,
    },
    "BILL_WENDY_CARROTS.TKA": {
        "detector": 3,
        "sample_mass_g": 38.6,
        "expected_live_time_s": 574753,
        "filled_height_cm": 0.91,
    },
    "CABBAGE.TKA": {
        "detector": 3,
        "sample_mass_g": 168.4,
        "expected_live_time_s": 13172,
        "filled_height_cm": 3.9,
    },
    "ESSEX LEEK.TKA": {
        "detector": 3,
        "sample_mass_g": 284.7,
        "expected_live_time_s": 574789,
        "filled_height_cm": 3.9,
    },
    "KENDAL GARLIC.TKA": {
        "detector": 3,
        "sample_mass_g": 90.8,
        "expected_live_time_s": 12872,
        "filled_height_cm": 3.9,
    },
    "Parsnip .TKA": {
        "detector": 3,
        "sample_mass_g": 173.8,
        "expected_live_time_s": 574793,
        "filled_height_cm": 1.2,
    },
}


BACKGROUND_INFO = {
    1: {
        "filename": "Background 1 WEEK.TKA",
        "expected_live_time_s": 575010,
    },
    2: {
        "filename": "HPGe2weekbackground.TKA",
        "expected_live_time_s": 575551,
    },
    3: {
        "filename": "WeekBackgroundDET3.TKA",
        "expected_live_time_s": 574748,
    },
}


# Detector calibration and efficiency fits from the project. Efficiency is
# stored as a fraction. These parameter errors do not cover every systematic.
DETECTOR_INFO = {
    1: {
        "slope": 0.2262,
        "offset": -122.31,
        "eff_coefficient": 36.5,
        "eff_power": -0.82,
        "calibration_uncertainty": 0.127,
        "eff_coefficient_uncertainty": 3.57,
        "eff_power_uncertainty": 0.015,
        "efficiency_unit": "fraction",
    },
    2: {
        "slope": 0.2762,
        "offset": 0.93,
        "eff_coefficient": 17.8,
        "eff_power": -0.85,
        "calibration_uncertainty": 0.125,
        "eff_coefficient_uncertainty": 2.93,
        "eff_power_uncertainty": 0.025,
        "efficiency_unit": "fraction",
    },
    3: {
        "slope": 0.2674,
        "offset": -146.03,
        "eff_coefficient": 15.4,
        "eff_power": -0.74,
        "calibration_uncertainty": 0.331,
        "eff_coefficient_uncertainty": 2.23,
        "eff_power_uncertainty": 0.022,
        "efficiency_unit": "fraction",
    },
}
