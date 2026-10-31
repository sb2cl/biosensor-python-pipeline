"""Mechanistic FdeR model families. See MODEL_LINEAGE.md before citing any result."""
from .minimal_fder import MinimalParameters, minimal_rhs, simulate_minimal_dose_response
from .final_fder import FinalModelData, load_final_model, simulate_final_dose_response
from .historical_2025 import HistoricalParameters, historical_rhs
