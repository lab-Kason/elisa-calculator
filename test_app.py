"""Headless UI smoke test for app.py (Streamlit AppTest). Run: python test_app.py"""

from __future__ import annotations

import pandas as pd
from streamlit.testing.v1 import AppTest

STD_COLS = ["Standard", "Concentration", "OD450", "OD600"]
SMP_COLS = ["Sample", "Dilution factor", "OD450", "OD600"]


def main() -> None:
    at = AppTest.from_file("app.py", default_timeout=120).run()
    assert len(at.exception) == 0, at.exception
    print("empty run OK")
    print("info:", [i.value[:70] for i in at.info])

    # Load the reference plate through session state (same path as a CSV upload)
    at.session_state["std_data"] = pd.read_csv("sample_data/example_A_clean.csv")[STD_COLS]
    at.session_state["smp_data"] = pd.read_csv("sample_data/example_A_samples_clean.csv")[SMP_COLS]
    at.session_state["std_ver"] = 1
    at.session_state["smp_ver"] = 1
    at.run()
    assert len(at.exception) == 0, at.exception
    print("reference-plate run OK")
    print("metrics:", [(m.label, m.value) for m in at.metric])
    print("dataframes:", len(at.dataframe))
    for frame in at.dataframe:
        print(frame.value.to_string())

    values = {m.label: m.value for m in at.metric}
    assert values.get("A · 零响应 (min)", "").startswith("0.08"), values
    assert values.get("D · 最大响应 (max)"), values
    results = at.dataframe[-1].value
    print("result table columns:", list(results.columns))
    assert "Reported conc" in " ".join(results.columns), results.columns
    reported = pd.to_numeric(results[results["Sample"] == "Wastewater original"]["Reported conc (pg/mL)"])
    print("original reported:", reported.tolist())

    # fresh session -> empty state (must not crash)
    at2 = AppTest.from_file("app.py", default_timeout=120).run()
    assert len(at2.exception) == 0, at2.exception
    print("cleared run OK")
    print("\nAppTest checks passed.")


if __name__ == "__main__":
    main()
