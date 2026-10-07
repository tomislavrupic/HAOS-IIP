# Semantic reporting repair

The first report mapped a failed per-level recovery-rate gate to `FIXED_BUDGET_NOT_SCALE_STABLE` even though the independently frozen scale-stability gate passed. The terminal mapping was corrected to `PARTIAL_RECOVERY_ONLY`. No raw row, seed, metric, threshold, gate, representation, decoder, control, or claim ceiling changed, and no final simulation was rerun.
