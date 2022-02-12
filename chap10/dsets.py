from collections import namedtuple
import csv
import functools
import glob
import os

import SimpleITK as sitk
import numpy as np


CandidateInfoTuple = namedtuple(
    "CandidateInfoTuple",
    "isNodule_bool, diameter_mm, series_uid, center_xyz",
)

@functools.lru_cache(1) # インメモリでキャッシングを行う
def getCandidateInfoList(requireOnDesk_bool=True):
    # 
    #
    #
    mhd_list = glob.glob("luna/subset*/*.mhd")
    presentOnDisk_set = {os.path.split(p)[-1][:-4] for p in mhd_list}

    diameter_dict = dict()
    with open("luna/annotation.csv", "r") as f:
        for row in list(csv.reader(f))[1:]:
            series_uid = row[0]
            annotationCenter_xyz = tuple([float(x) for x in row[1:4]])
            annotationDiameter_mm = float(row[4])

            diameter_dict.setdefault(series_uid, []).append(
                (annotationCenter_xyz, annotationDiameter_mm)
            )
    
    candidateInfo_list = list()
    with open("luna/candidates.csv", "r") as f:
        for row in list(csv.reader(f))[1:]:
            series_uid = row[0]

            if series_uid not in presentOnDisk_set and requireOnDesk_bool:
                continue

            isNodule_bool = bool(int(row[4]))
            candidateCenter_xyz = tuple([float(x) for x in row[1:4]])

            candidateDiameter_mm = 0.0
            for annotatin_tup in diameter_dict.get(series_uid, []):
                annotationCenter_xyz, annotationDiameter_mm = annotatin_tup
                for i in range(3):
                    delta_mm = abs(candidateCenter_xyz[i] - annotationCenter_xyz[i])
                    if delta_mm > annotationDiameter_mm / 4:
                        break # 中心座標のずれが大きい場合はbreak
                    else:
                        candidateDiameter_mm = annotationDiameter_mm
                        break
            
            candidateInfo_list.append(CandidateInfoTuple(
                isNodule_bool,
                candidateDiameter_mm,
                series_uid,
                candidateCenter_xyz,
            ))
    
    candidateInfo_list.sort(reverse=True)
    return candidateInfo_list

class Ct:
    def __init__(self, series_uid):
        mhd_path = glob.glob(
            "luna/subset*/{}.mhd".format(series_uid)
        )[0]

        ct_mhd = sitk.ReadImage(mhd_path)
        ct_a = np.array(sitk.GetArrayFromImage(ct_mhd), dtype=np.float32)

        #
        #
        #
        ct_a.clip(-1000, 1000, ct_a)

        self.seried_uid = series_uid
        self.hu_a = ct_a