"""Unit tests for the [ANFREQ] placeholder that replaced hard-coded "! AnFreq".

The refactor's whole claim is that it changed *nothing* by default: every mode
renders the same AnFreq state its template used to hard-code. These tests pin
that claim to the pre-refactor truth, recorded here as an explicit table, so a
future edit to a template or to MODE_ANFREQ_DEFAULT cannot quietly flip a mode's
default on or off.
"""

from __future__ import annotations

import pytest

from xas_pipeline import resources, templates
from xas_pipeline.stages import orca_prep

# What each template hard-coded before [ANFREQ] existed (git 6b970e2). Written out
# rather than read from MODE_ANFREQ_DEFAULT: a test that compares the table to
# itself would pass no matter what the table said.
ANFREQ_BEFORE_REFACTOR = {
    "caopt-anfreq": True,
    "quick": False,
    "quick-ca-fixed": False,
    "hopt-anfreq": True,
    "carved-anfreq": True,
    "no-constraints": True,
    "backbone": True,
    "xtb-free": True,
    "xtb-constrained": False,
    "carved-spring": False,
    "hopt-spring": False,
}


def test_every_registered_mode_has_a_default():
    assert set(orca_prep.MODE_ANFREQ_DEFAULT) == set(orca_prep.TEMPLATE_FILE_BY_MODE)


def test_defaults_match_what_the_templates_used_to_hard_code():
    assert orca_prep.MODE_ANFREQ_DEFAULT == ANFREQ_BEFORE_REFACTOR


@pytest.mark.parametrize("mode", sorted(orca_prep.TEMPLATE_FILE_BY_MODE))
def test_template_carries_both_placeholders(mode):
    """A template missing [ANFREQ] would silently ignore --anfreq/--no-anfreq, and
    one missing [GEOMETRY] would ignore a staged run's starting geometry."""
    text = (resources.template_root() / orca_prep.TEMPLATE_FILE_BY_MODE[mode]).read_text()
    assert "[ANFREQ]" in text, f"{mode}: template has no [ANFREQ] placeholder"
    assert "[GEOMETRY]" in text, f"{mode}: template has no [GEOMETRY] placeholder"
    # The keyword must not ALSO be hard-coded, or --no-anfreq could not remove it.
    assert "! AnFreq" not in text, f"{mode}: template still hard-codes ! AnFreq"


@pytest.mark.parametrize("mode", sorted(orca_prep.TEMPLATE_FILE_BY_MODE))
def test_rendering_with_the_default_reproduces_the_old_keyword_state(mode):
    text = (resources.template_root() / orca_prep.TEMPLATE_FILE_BY_MODE[mode]).read_text()
    enabled = orca_prep.anfreq_default_for(mode)
    rendered = templates.fill(text, {"ANFREQ": orca_prep.anfreq_directive(enabled)})
    directives = [ln for ln in rendered.splitlines() if ln.strip() == "! AnFreq"]
    assert directives == (["! AnFreq"] if ANFREQ_BEFORE_REFACTOR[mode] else [])


@pytest.mark.parametrize("mode", sorted(orca_prep.TEMPLATE_FILE_BY_MODE))
def test_the_toggle_can_move_every_mode_both_ways(mode):
    text = (resources.template_root() / orca_prep.TEMPLATE_FILE_BY_MODE[mode]).read_text()
    on = templates.fill(text, {"ANFREQ": orca_prep.anfreq_directive(True)})
    off = templates.fill(text, {"ANFREQ": orca_prep.anfreq_directive(False)})
    assert "\n! AnFreq\n" in on
    assert "! AnFreq" not in off


def test_plan_stages_single_stage_uses_the_mode_default():
    assert orca_prep.plan_stages("caopt-anfreq", [], None) == [("caopt-anfreq", True)]
    assert orca_prep.plan_stages("quick-ca-fixed", [], None) == [("quick-ca-fixed", False)]


def test_plan_stages_override_applies_to_the_final_stage_only():
    stages = orca_prep.plan_stages("caopt-anfreq", ["quick-ca-fixed"], False)
    assert stages == [("quick-ca-fixed", False), ("caopt-anfreq", False)]

    stages = orca_prep.plan_stages("quick-ca-fixed", ["quick"], True)
    assert stages == [("quick", False), ("quick-ca-fixed", True)]


def test_pre_stages_never_run_anfreq_even_when_their_mode_defaults_on():
    """An intermediate geometry is about to move again, so its Hessian is waste."""
    stages = orca_prep.plan_stages("caopt-anfreq", ["carved-anfreq", "hopt-anfreq"], None)
    assert [af for _, af in stages] == [False, False, True]
