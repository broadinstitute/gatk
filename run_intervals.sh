./scripts/tf.sh -t \
    /home/${USER}/repos/ml/ml4cvd/recipes.py \
    --mode train \
    --tensors /home/${USER}/partners_ecg/hd5_copies \
    --input_tensors partners_ecg_voltage \
    --output_tensors \
                        partners_ecg_rate \
                        partners_ecg_qrs \
                        partners_ecg_pr \
                        partners_ecg_qt \
                        partners_ecg_qtc \
    --inspect_model \
    --epochs 300 \
    --batch_size 256 \
    --training_steps 256 \
    --validation_steps 32 \
    --test_steps 16 \
    --patience 10 \
    --test_modulo 0 \
    --learning_rate 0.0002 \
    --output_folder "/home/${USER}/Dropbox\ \(Partners\ HealthCare\)/partners_ecg/ml4cvd_results/" \
    --id intervals_all_not_norm
