## pfamily_list validation shared with glmbayesCore::multi_rlmb (internal API).
## Resolved at load time so block_lmb / block_glmb stay aligned with multi-response
## pfamily conventions without duplicating validation logic.

.mrglmb_normalize_pfamily_lists <- getFromNamespace(
  ".mrglmb_normalize_pfamily_lists", "glmbayesCore"
)
.validate_pfamily_for_rlmb <- getFromNamespace(
  ".validate_pfamily_for_rlmb", "glmbayesCore"
)
