from pathlib import Path
import pandas as pd
rows=[]
for out in sorted(Path('.').glob('*.out')):
    text=out.read_text(errors='ignore')
    ok='ORCA TERMINATED NORMALLY' in text
    in_ir=False
    for line in text.splitlines():
        if 'IR SPECTRUM' in line:
            in_ir=True
            continue
        if in_ir and ('The first frequency' in line or 'THERMOCHEMISTRY' in line or 'NORMAL MODES' in line):
            break
        parts=line.split()
        if len(parts)>=3:
            try:
                mode=int(parts[0].replace(':',''))
                freq=float(parts[1])
                nums=[]
                for p in parts[1:]:
                    try: nums.append(float(p))
                    except Exception: pass
                intensity=nums[-1] if nums else None
                if 450 <= freq <= 1200:
                    rows.append({'file':out.name,'terminated_normally':ok,'mode':mode,'frequency_cm-1':freq,'IR_intensity_or_last_numeric_column':intensity})
            except Exception:
                pass
pd.DataFrame(rows).to_csv('orca_ir_modes_450_1200.csv',index=False)
print(pd.DataFrame(rows))
print('Saved: orca_ir_modes_450_1200.csv')
