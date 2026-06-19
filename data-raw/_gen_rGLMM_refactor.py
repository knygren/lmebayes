from pathlib import Path
ROOT = Path(r"c:\Rpackages\lmebayes\R")

RGLMM_TEMP_PATH = Path(r"c:\Rpackages\lmebayes\data-raw\_rGLMM_temp_body.R")
RGLMERB_PATH = Path(r"c:\Rpackages\lmebayes\data-raw\_rglmerb_body.R")

ROOT.joinpath("rGLMM_temp.R").write_text(RGLMM_TEMP_PATH.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
ROOT.joinpath("rglmerb.R").write_text(RGLMERB_PATH.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
print("done")
