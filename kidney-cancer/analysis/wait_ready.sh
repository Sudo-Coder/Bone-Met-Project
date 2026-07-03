#!/bin/bash
cd /autofs/projects-t3/hussain/scProj
for i in $(seq 1 90); do
  env_ok=0; dl_ok=0
  envs/mechanism_env/bin/python -c "import pyscenic" 2>/dev/null && env_ok=1
  grep -q "DONE downloads" outputs/logs/download_refs.log 2>/dev/null && dl_ok=1
  if [ $env_ok -eq 1 ] && [ $dl_ok -eq 1 ]; then
     echo "READY env=$env_ok dl=$dl_ok after $((i*20))s"; break
  fi
  sleep 20
done
echo "=== env ==="; envs/mechanism_env/bin/python -c "import pyscenic,decoupler; print('pyscenic',pyscenic.__version__,'decoupler',decoupler.__version__)" 2>&1 | tail -1
echo "=== refs ==="; ls -lh resources/cistarget resources/nichenet 2>/dev/null | grep -vE '^total|^d'
grep -E "Successfully installed" outputs/logs/build_env.log >/dev/null && echo "pip: installed" || echo "pip: NOT done"
