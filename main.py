import csv
import pandas as pd

HMEC_dict = {}

interactions_forGSM = {}
Non_overlaped_interactions={}
# overlappedPromoter = {}
overlappedGSM_with_promoterInteraction ={}
overlappedGSM_with_NonOverlapedInteraction=[]
dataframe_GSM=pd.DataFrame()

addresOHmecFile = "h.csv"
df_start_hmec = pd.read_csv(addresOHmecFile, usecols=['start'], low_memory=True)
df_end_hmec = pd.read_csv(addresOHmecFile, usecols=['end'], low_memory=True)
df_chr_hmec = pd.read_csv(addresOHmecFile, usecols=["chr"], low_memory=True)
df_start_of_interaction = pd.read_csv(addresOHmecFile, usecols=["start_interaction"], low_memory=True)
df_end_of_interaction = pd.read_csv(addresOHmecFile, usecols=["end_interaction"], low_memory=True)
df_chr_startHmec_startInteraction = pd.read_csv(addresOHmecFile, usecols=["chr", "start", "start_interaction"],low_memory=True)

def add_Required_HmecData_toDictionary(dataFrame1, dic, dataframe2, dataframe3):
    for i in range(len(dataFrame1.index)):
        key = int(dataFrame1.iloc[i].to_string(index=False))
        if key not in dic.keys():
            start_end_of_intraction = ( int(dataframe2.iloc[i].to_string(index=False)),
                                       int(dataframe3.iloc[i].to_string(index=False)) )

            dic[key] = [int(df_end_hmec.iloc[i].to_string(index=False))
                , int(df_chr_hmec.iloc[i].to_string(index=False)), start_end_of_intraction]
            del start_end_of_intraction

        else:
            key = int(dataFrame1.iloc[i].to_string(index=False))
            start_end_of_intraction = ( int(dataframe2.iloc[i].to_string(index=False)),
                                       int(dataframe3.iloc[i].to_string(index=False)) )
            dic[key].append(start_end_of_intraction)
            del start_end_of_intraction

add_Required_HmecData_toDictionary(df_start_hmec, HMEC_dict, df_start_of_interaction, df_end_of_interaction)
print(HMEC_dict)



# def hash_start_end(num):
#     return int(int(num) // 5000) * 5000
#
# numberOf_allPromoters = 0
#
# with open('p.csv', 'r') as Promfile:
#     csv_reader = csv.reader(Promfile)
#     dataframe_promoter = pd.DataFrame(csv_reader)
#     index = 0
#     numberOf_allPromoters = (len(dataframe_promoter.index)) - 1
#
#     for i in range(1, len(dataframe_promoter.index) ):
#         begin_overflow = dataframe_promoter.iloc[[i]][1].to_string(index=False)
#         end_overflow = dataframe_promoter.iloc[[i]][2].to_string(index=False)
#         start = hash_start_end(int(begin_overflow))  # hash shode
#         end = hash_start_end(int(end_overflow))  # hash shode
#
#         if start == end:
#             if start in HMEC_dict:
#                 HMEC_dict[start].append("-")
#                 HMEC_dict[start].append([f"promoter_{str(index)}", (begin_overflow), end_overflow])
#
#
#         else:
#             if start in HMEC_dict:
#                 if int(start) + 5000 - int(begin_overflow) >= 10:
#                     HMEC_dict[start].append("-")
#                     HMEC_dict[start].append([f"promoter_{str(index)}", begin_overflow, end_overflow])
#
#             if end in HMEC_dict:
#                 if int(end_overflow) - int(end) >= 10:
#                     HMEC_dict[end].append("-")
#                     HMEC_dict[end].append([f"promoter_{str(index)}", begin_overflow, end_overflow])
#
#
#         index = index + 1
#
#
# # print(HMEC_dict)  {710000: [715000, 9, (224180000, 224185000), (3445, 56), '-', ['promoter_91', '712664', '714664'],    760000: [(9999, 715000, 760000)]}
#
# allData = []
#
# def find_overlappedPromoter_with_Hmecs(HmecDataframe, HmecDict, alldataArray ):
#     for i in range(len(HmecDataframe.index)):  # promoter
#         if int(HmecDataframe.iloc[i].to_string(index=False)) in HMEC_dict.keys():
#             key = int(HmecDataframe.iloc[i].to_string(index=False))
#             value = HmecDict[key]
#             if '-' in value:    #فقظ اونایی رو که پروموتر دارند اد کن
#                 dual = (key, value)
#                 if dual not in alldataArray:
#                     alldataArray.append(dual)
#                 del dual
#
# find_overlappedPromoter_with_Hmecs(df_start_hmec, HMEC_dict, allData )
#
# headerList = ['start of HIC','[ end of HIC , chr_hmec , ( start of interactions , end of interactions )  , overlap with = [promoterID , start of promoter , end of promoter] ']
# if allData != []:
#     with open('comparedFilePromoter.csv', 'w') as csv_file:
#         writer = csv.writer(csv_file)
#         writer.writerow(headerList)
#         for i in allData:
#             writer.writerow(i)
#
# # ----------------------------------------------------------------------------------------------------------------------
#
# numberOf_matchedInteraction = 0
# def find_overlappedInteraction():
#     for i in allData:
#         start_hmec = i[0]
#         charr = i[1][1]
#         for j in i[1]:
#             if type(j) == tuple:
#                 startOF_interaction = j[0]
#                 interactions_forGSM[startOF_interaction] = [((charr, start_hmec, startOF_interaction))]
#             # interactions_forGSM اینترکشن های راه پیدا کرده به مرحله بعد
#
# find_overlappedInteraction()
#
#
# dataframe_GSM=pd.DataFrame()
# numberOf_allGSM = 0
# with open('g.csv', 'r') as GSMfile:
#     csv_reader = csv.reader(GSMfile)
#     dataframe_GSM = pd.DataFrame(csv_reader)
#     numberOf_allGSM = (len(dataframe_GSM.index)) - 1
#
#
#     for i in range(1, (len(dataframe_GSM.index))):
#         begin_of_GSM = dataframe_GSM.iloc[[i]][1].to_string(index=False)
#         end_of_GSM = dataframe_GSM.iloc[[i]][2].to_string(index=False)
#         start_hashed_GSM = hash_start_end(int(begin_of_GSM))  # hash shode
#         end_hashed_GSM = hash_start_end(int(end_of_GSM))  # hash shode
#
#         if start_hashed_GSM == end_hashed_GSM:
#             if start_hashed_GSM in interactions_forGSM:
#                 interactions_forGSM[start_hashed_GSM].append([begin_of_GSM, end_of_GSM])
#                 overlappedGSM_with_promoterInteraction[begin_of_GSM]=0  #فقط برای اینکه بفهمم چه جی اس ام هایی اورلپ دارن
#
#         else:
#             if start_hashed_GSM in interactions_forGSM:
#                 if int(start_hashed_GSM) + 5000 - int(begin_of_GSM) >= 10:
#                     interactions_forGSM[start_hashed_GSM].append([begin_of_GSM, end_of_GSM])
#                     overlappedGSM_with_promoterInteraction[begin_of_GSM]=0
#
#             if end_hashed_GSM in interactions_forGSM:
#                 if int(end_of_GSM) - int(end_hashed_GSM) >= 10:
#                     interactions_forGSM[end_hashed_GSM].append([begin_of_GSM, end_of_GSM])  #برای اینترکشن های پیدا شده میایم ببینیم چه جی اس ام هایی در این بازه ها اورلب دارند
#                     overlappedGSM_with_promoterInteraction[begin_of_GSM]=0  #اون جی اس ام هایی که اولپ داره  رو نگه میداره
#
#
#
#
# # for removing interactions that didnt overlaped
# def remove_overlappedInteraction_fromDataFrame(dict, dataframe):
#     for k, v in dict.items():
#         chr, hmec, interaction = v[0]
#         index_names = dataframe[ (dataframe['chr'] == chr) & (dataframe['start'] == hmec) ].index
#         dataframe.drop(index_names, inplace=True)
#     return (len(dataframe))
#
# numberof_nonOverlappedInteraction=remove_overlappedInteraction_fromDataFrame(interactions_forGSM, df_chr_startHmec_startInteraction)
#
#
# def find_nonOverlaped_interactions(dataframe, Non_overlaped_interactionsDict):
#     for i in range(len(dataframe.index)):
#         key = (dataframe['start_interaction'].iloc[i])
#         chr_startHmec=(dataframe['start'].iloc[i] , dataframe['chr'].iloc[i]  )
#         if key not in Non_overlaped_interactionsDict.keys():
#             Non_overlaped_interactionsDict[key] = [chr_startHmec]
#             del chr_startHmec
#         else:
#             key = (dataframe['start_interaction'].iloc[i])
#             chr_startHmec=(dataframe['start'].iloc[i] , dataframe['chr'].iloc[i]  )
#             Non_overlaped_interactionsDict[key].append(chr_startHmec)
#             del chr_startHmec
#
# find_nonOverlaped_interactions(df_chr_startHmec_startInteraction, Non_overlaped_interactions)
#
#
#
#
#
#
# oneof_their_heads_contains_GSM=[]
# for i in range(1, (len(dataframe_GSM.index))):
#     begin_of_GSM = dataframe_GSM.iloc[[i]][1].to_string(index=False)
#     end_of_GSM = dataframe_GSM.iloc[[i]][2].to_string(index=False)
#     start_hashed_GSM = hash_start_end(int(begin_of_GSM))  # hash shode
#     end_hashed_GSM = hash_start_end(int(end_of_GSM))  # hash shode
#
#     if start_hashed_GSM == end_hashed_GSM:
#         if start_hashed_GSM in Non_overlaped_interactions:
#             oneof_their_heads_contains_GSM.append([begin_of_GSM, end_of_GSM])
#
#
#     else:
#         if start_hashed_GSM in Non_overlaped_interactions:
#             if int(start_hashed_GSM) + 5000 - int(begin_of_GSM) >= 10:
#                 oneof_their_heads_contains_GSM.append([begin_of_GSM, end_of_GSM])
#
#         if end_hashed_GSM in Non_overlaped_interactions:
#             if int(end_of_GSM) - int(end_hashed_GSM) >= 10:
#                 oneof_their_heads_contains_GSM.append([begin_of_GSM, end_of_GSM])
#
#
#
#
#
# df = pd.DataFrame.from_dict(interactions_forGSM, orient="index")
# df.to_csv("comparedFileGSM.csv")
#
#
#
# def calculate_division(first, denominator):
#     answer=0
#     if (denominator) !=0:
#         answer= float( (first)/(denominator) )*100
#         return answer
#     else:
#         return answer
#
#
# def write_information():
#
#     lenOf_Hmec = len(df_start_hmec.index)
#
#     dict = {"Number of all interactions ": lenOf_Hmec ,
#             "number of all promoters ": numberOf_allPromoters,
#             "Number of all GSMs ": numberOf_allGSM,
#
#             'Number of overlaped interaction (just contain promoter) ': (len(interactions_forGSM)),#(یک سر پروموتر دار) تعداد تمام اینترکشن هایی که به مرحله بعد آمدند
#
#             "Number of GSMs that had overlap with interactions (contain promoter and GSM) ":len(overlappedGSM_with_promoterInteraction) , #   یه سر پروموتر یه سر جی اس ام
#
#             "Number of nonOverlap interactions ":numberof_nonOverlappedInteraction ,#تعداد اینترکشن هایی که به مرحله بعد نیومدن (فاقد پروموتر)
#
#             "Number of Non-overlapped interaction that some GSMs are in their bin ":len(oneof_their_heads_contains_GSM),  #تعداد اینترکشن هایی که به مرحله بعد نیومدن اما این سری از اینها با جی اس ام ها اورلپ دارند
#
#     }
#     print(Non_overlaped_interactions)
#     print(oneof_their_heads_contains_GSM)
#
#
#
#     with open("information.csv", 'w') as csv_file:
#         writer = csv.writer(csv_file)
#         for key, value in dict.items():
#             writer.writerow([key, value])
#
# write_information()













# import pandas as pd
# import matplotlib.pyplot as plt
#
# file_path = "/information.csv"
#
# one_side_promoter_overlapped_precent = []
# no_side_promoter_overlapped_precent = []
#
# for i in range(1,12):
#
#     data = pd.read_csv(str(i)+file_path)
#     promoter_interaction = int(list(data.iloc[1])[1])
#     histon_promoter_interaction = int(list(data.iloc[6])[1])
#     nonpromoter_interaction = int(list(data.iloc[9])[1])
#     nonpromoter_histon_interaction = int(list(data.iloc[10])[1])
#
#     percentage_promoter = histon_promoter_interaction / promoter_interaction
#     percentage_nonpromoter = nonpromoter_histon_interaction / nonpromoter_interaction
#     one_side_promoter_overlapped_precent.append(percentage_promoter)
#     no_side_promoter_overlapped_precent.append(percentage_nonpromoter)
#
# result_df = pd.DataFrame({
#     'labels': [i for i in range(1,12)],
#     'promo_overlapped': one_side_promoter_overlapped_precent,
#     'non_promo_overlapped': no_side_promoter_overlapped_precent
# })
#
# print(result_df.head())
# result_df.plot("labels", ['promo_overlapped', 'non_promo_overlapped'], kind='bar')
# plt.show()
#


