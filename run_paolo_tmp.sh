bash ./scripts/tf_paolo.sh -t \
    /home/${USER}/repos/ml/ml4cvd/recipes.py \
    --mode plot_histograms \
    --tensors /home/${USER}/partners_ecg/hd5/2019-01 \
    --input_tensors partners_ecg_voltage \
    --output_tensors \
                     partners_ecg_qt \
    --inspect_model \
    --epochs 5 \
    --batch_size 128 \
    --training_steps 8192 \
    --validation_steps 128 \
    --test_steps 64 \
    --patience 10 \
    --test_modulo 0 \
    --learning_rate 0.00002 \
    --output_folder "/home/${USER}/paolo_test_intervals/" \
    --id qt_test
