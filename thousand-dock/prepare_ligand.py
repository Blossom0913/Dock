import os
import subprocess
import sys
import time

# 配置路径
BASE_DIR = r"D:\2025Spring\protein-predic\thousand-dock"

OBABEL_PATH = r"D:\2025Spring\protein-predic\OpenBabel-3.1.1\obabel.exe"  # 修改为你的实际路径

LIGANDS_DIR = os.path.join(BASE_DIR, "ligands")
PREPARED_DIR = os.path.join(BASE_DIR, "prepared")
LOGS_DIR = os.path.join(BASE_DIR, "logs")

# 确保目录存在
os.makedirs(PREPARED_DIR, exist_ok=True)
os.makedirs(LOGS_DIR, exist_ok=True)

# 设置MGLTools路径（根据实际安装路径修改）
MGLTOOLS_PATH = r"D:\2025Spring\mgltools"
PREPARE_LIGAND = os.path.join(
    MGLTOOLS_PATH, "MGLToolsPckgs", "AutoDockTools", "Utilities24", "prepare_ligand4.py"
)
PYTHON_EXE = os.path.join(MGLTOOLS_PATH, "python.exe")

# 日志文件
LOG_FILE = os.path.join(LOGS_DIR, "ligand_preparation.log")
ERROR_LOG = os.path.join(LOGS_DIR, "preparation_errors.log")

def log(message, file=LOG_FILE):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(file, "a") as f:
        f.write(f"[{timestamp}] {message}\n")
    print(message)

# 清空日志文件
open(LOG_FILE, "w").close()
open(ERROR_LOG, "w").close()

log("Starting ligand preparation process")
log(f"Ligands directory: {LIGANDS_DIR}")
log(f"Prepared directory: {PREPARED_DIR}")

# 计数器
success_count = 0
fail_count = 0

# 处理每个配体文件
for filename in os.listdir(LIGANDS_DIR):
    if filename.lower().endswith(".pdb"):
        ligand_name = os.path.splitext(filename)[0]
        input_pdb = os.path.join(LIGANDS_DIR, filename)
        temp_pdb = os.path.join(PREPARED_DIR, f"{ligand_name}_temp.pdb")
        output_pdbqt = os.path.join(PREPARED_DIR, f"{ligand_name}.pdbqt")
        
        log(f"Processing: {ligand_name}")
        
        try:
            # 步骤1: 使用OpenBabel预处理（加氢、修复结构）
            obabel_cmd = f'"{OBABEL_PATH}" "{input_pdb}" -O "{temp_pdb}" -h'
            subprocess.run(obabel_cmd, shell=True, check=True, capture_output=True)
            
            # 替代 prepare_ligand4.py 的步骤
            obabel_cmd2 = f'"{OBABEL_PATH}" "{temp_pdb}" -O "{output_pdbqt}"'
            subprocess.run(obabel_cmd2, shell=True, check=True, capture_output=True)
            result = subprocess.run(
                obabel_cmd2, shell=True, check=True, capture_output=True, text=True
            )
            
            # 检查输出文件
            if os.path.exists(output_pdbqt) and os.path.getsize(output_pdbqt) > 0:
                success_count += 1
                log(f"Success: {ligand_name} converted to PDBQT")
                # 删除临时文件
                os.remove(temp_pdb)
            else:
                raise Exception("PDBQT file not created or empty")
                
        except subprocess.CalledProcessError as e:
            error_msg = f"Error processing {ligand_name}: {e.stderr}"
            log(error_msg, ERROR_LOG)
            log(f"Failed: {ligand_name} - {str(e)}")
            fail_count += 1
        except Exception as e:
            error_msg = f"General error with {ligand_name}: {str(e)}"
            log(error_msg, ERROR_LOG)
            log(f"Failed: {ligand_name} - {str(e)}")
            fail_count += 1

# 生成报告
log("\n===== PREPARATION REPORT =====")
log(f"Total ligands processed: {success_count + fail_count}")
log(f"Successfully prepared: {success_count}")
log(f"Failed: {fail_count}")
log("Process completed")

if fail_count > 0:
    log(f"Check error log for details: {ERROR_LOG}")