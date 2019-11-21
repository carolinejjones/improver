# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# (C) British Crown Copyright 2017-2019 Met Office.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
"""Tests for the wxcode CLI"""

import argparse
import pytest

from improver.cli import wxcode
from improver.tests import acceptance as acc


@pytest.mark.acc
@acc.skip_if_kgo_missing
def test_basic(tmp_path):
    """Test basic wxcode processing"""
    kgo_dir = acc.kgo_root() / "wxcode/basic"
    kgo_path = kgo_dir / "kgo.nc"
    params = ["lightning_flashes_per_unit_area_in_vicinity_above",
              "cloud_area_fraction_above",
              "rainfall_rate_above",
              "rainfall_rate_in_vicinity_above",
              "lwe_snowfall_rate_above",
              "lwe_snowfall_rate_in_vicinity_above",
              "visibility_in_air_below",
              "low_type_cloud_area_fraction_above"]
    param_paths = [str(kgo_dir / f"probability_of_{p}_threshold.nc")
                   for p in params]
    output_path = tmp_path / "output.nc"
    args = [*param_paths, str(output_path)]
    wxcode.main(args)
    acc.compare(output_path, kgo_path)


@pytest.mark.acc
@acc.skip_if_kgo_missing
def test_native_units(tmp_path):
    """Test wxcode processing with native units"""
    kgo_dir = acc.kgo_root() / "wxcode/basic"
    input_dir = acc.kgo_root() / "wxcode/native_units"
    kgo_path = kgo_dir / "kgo.nc"
    params = ["cloud_area_fraction_above",
              "lightning_flashes_per_unit_area_in_vicinity_above",
              "rainfall_rate_above",
              "rainfall_rate_in_vicinity_above",
              "lwe_snowfall_rate_above",
              "lwe_snowfall_rate_in_vicinity_above",
              "visibility_in_air_below",
              "low_type_cloud_area_fraction_above"]
    param_paths = [str(input_dir / f"probability_of_{p}_threshold.nc")
                   for p in params]
    output_path = tmp_path / "output.nc"
    args = [*param_paths, str(output_path)]
    wxcode.main(args)
    acc.compare(output_path, kgo_path)


@pytest.mark.acc
@acc.skip_if_kgo_missing
def test_global(tmp_path):
    """Test global wxcode processing"""
    kgo_dir = acc.kgo_root() / "wxcode/global"
    kgo_path = kgo_dir / "kgo.nc"
    params = ["rainfall_rate_above",
              "snowfall_rate_above",
              "visibility_at_screen_level_below",
              "cloud_area_fraction_above",
              "low_type_cloud_area_fraction_above"]
    param_paths = [str(kgo_dir / f"probability_of_{p}_threshold.nc")
                   for p in params]
    output_path = tmp_path / "output.nc"
    args = ["--wxtree=global", *param_paths, str(output_path)]
    wxcode.main(args)
    acc.compare(output_path, kgo_path)


@pytest.mark.acc
@acc.skip_if_kgo_missing
def test_insufficent_files(tmp_path):
    """Test wxcode processing with insufficent files"""
    kgo_dir = acc.kgo_root() / "wxcode/global"
    params = ["rainfall_rate_above",
              "snowfall_rate_above",
              "cloud_area_fraction_above",
              "low_type_cloud_area_fraction_above"]
    param_paths = [str(kgo_dir / f"probability_of_{p}_threshold.nc")
                   for p in params]
    output_path = tmp_path / "output.nc"
    args = ["--wxtree=global", *param_paths, str(output_path)]
    with pytest.raises(OSError):
        wxcode.main(args)


@pytest.mark.acc
@acc.skip_if_kgo_missing
def test_no_lightning(tmp_path):
    """Test wxcode processing with no lightning"""
    kgo_dir = acc.kgo_root() / "wxcode/basic"
    kgo_path = kgo_dir / "kgo_no_lightning.nc"
    params = ["rainfall_rate_above",
              "rainfall_rate_in_vicinity_above",
              "lwe_snowfall_rate_above",
              "lwe_snowfall_rate_in_vicinity_above",
              "visibility_in_air_below",
              "cloud_area_fraction_above",
              "low_type_cloud_area_fraction_above"]
    param_paths = [str(kgo_dir / f"probability_of_{p}_threshold.nc")
                   for p in params]
    output_path = tmp_path / "output.nc"
    args = [*param_paths, str(output_path)]
    wxcode.main(args)
    acc.compare(output_path, kgo_path)
