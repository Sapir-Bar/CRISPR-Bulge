import os
os.environ['PYTHONHASHSEED'] = '0'

import pandas as pd
from OT_deep_score_src.evaluate_utilities import measure_clevage_acc
from OT_deep_score_src.general_utilities import Model_task, Model_type
from train_and_predict_scripts import predict_config

def ensemble_predict_wrapper():
    for ensemble_componet_i in range(5):
            predict_config.main(
                version="5_revision_normalized_ensemble_{}".format(ensemble_componet_i), setting_number=5,
                model_types=(Model_type.C_3,),
                model_tasks=(Model_task.REGRESSION_TASK,))
    
def measure_clevage_acc_wrapper(prediction_file_name):
    return measure_clevage_acc(
        prediction_file_name=prediction_file_name,
        models=["c_3-thReads-0-Regression-seq", ],
        only_bulges=True,
        only_mismatches=False,
        model_task=Model_task.REGRESSION_TASK
    )

def measure_clevage_acc_per_ensemble_wrapper(prediction_file_name):
    return measure_clevage_acc(
        prediction_file_name=prediction_file_name,
        models=[f"c_3-thReads-0-Regression-seq_ensemble_{i}" for i in range(5)],
        only_bulges=False,
        only_mismatches=False,
        model_task=Model_task.REGRESSION_TASK
    )

def concat_ensemble_results():
    pred_col = "c_3-thReads-0-Regression-seq"
    base_df = None

    for i in range(5):
        filename = (
            f"predictions_NewGUIDEseq_cleavage_v5_revision_normalized_ensemble_{i}_"
            f"exclude_NewGUIDEseq_continue_from_change_seq_folds_1.csv"
        )
        df = pd.read_csv(filename)

        if i == 0:
            base_df = df.copy()
            base_df.rename(
                columns={pred_col: f"{pred_col}_ensemble_{i}"},
                inplace=True
            )
        else:
            base_df[f"{pred_col}_ensemble_{i}"] = df[pred_col]

    output_file = "predictions_5_revision_normalized_ensemble_final.csv"
    base_df.to_csv(output_file, index=False)

    return output_file