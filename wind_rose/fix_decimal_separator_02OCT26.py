# fix_decimal_separator_02OCT26.py
# Repairs the Risopatron Campbell CR1000X weather-station export (sheet "Periodo 1", 22-31 January 2022). In that
# sheet, values written with three decimals ("2.367") were read with "." as a thousands separator and stored as
# integers (2367, number format "#,##0"). The logger writes four significant digits, so only values below 10 carry
# three decimals and those cells are divided by 1000; pressure is written at full precision, so affected pressure
# cells are divided by the power of 10 that returns them to 900-1100 hPa. Sheet "Periodo 2" is copied unchanged.
# Input  : DatosAll_CR1000X_RPT.xlsx (weather-station record, not included in this repository)
# Output : DatosAll_CR1000X_RPT_decimalfix_02OCT26.csv (one row per 5 min record)
# Run    : python3 fix_decimal_separator_02OCT26.py (from the wind_rose/ folder)

import math
import openpyxl
import pandas as pd

SRC = "DatosAll_CR1000X_RPT.xlsx"
OUT = "DatosAll_CR1000X_RPT_decimalfix_02OCT26.csv"
FP2 = {"WS_ms_S_WVT", "WindDir_D1_WVT", "AirTC_Avg", "RH_Avg", "BattV_Avg", "SlrW_Avg", "Raw_mV_Avg",
       "CS320_Temp_Avg", "CS320_X_Avg", "CS320_Y_Avg", "CS320_Z_Avg", "DewPtC_Avg"}

wb = openpyxl.load_workbook(SRC, read_only=True)
frames, n_fixed = [], {}
for sheet in ["Periodo 1", "Periodo 2"]:
    rows = list(wb[sheet].iter_rows(min_row=1))
    names = [c.value for c in rows[1]]
    recs = []
    for row in rows[4:]:
        rec = {}
        for name, c in zip(names, row):
            v = c.value
            if isinstance(v, (int, float)) and c.number_format == "#,##0" and name not in ("RECORD",):
                if name == "BP_mbar_Avg":
                    k = 0
                    while abs(v) / 10 ** k > 1100:
                        k += 1
                    v = v / 10 ** k
                elif name in FP2:
                    v = v / 1000
                n_fixed[(sheet, name)] = n_fixed.get((sheet, name), 0) + 1
            rec[name] = v
        rec["sheet"] = sheet
        recs.append(rec)
    frames.append(pd.DataFrame(recs))

d = pd.concat(frames, ignore_index=True)
d["TIMESTAMP"] = pd.to_datetime(d["TIMESTAMP"], errors="coerce")
d = d.dropna(subset=["TIMESTAMP"]).drop_duplicates("TIMESTAMP").sort_values("TIMESTAMP")
for c in names:
    if c != "TIMESTAMP":
        d[c] = pd.to_numeric(d[c], errors="coerce")

print("Cells repaired (sheet, column): count")
for k, v in sorted(n_fixed.items()):
    print(f"  {k}: {v}")
d.to_csv(OUT, index=False)
print("Wrote", OUT, len(d), "records")
