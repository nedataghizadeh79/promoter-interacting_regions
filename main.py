
import csv
import pandas as pd

HMEC_dict = {1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [],
             14: [], 15: [], 16: [], 17: [], 18: [], 19: [], 20: [], 21: [], 22: [], 23: [], 24: [], 25: []}

all_chars = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

interactions_forGSM = {}
Non_overlaped_interactions = {}
# overlappedPromoter = {}
overlappedGSM_with_promoterInteraction = {}
overlappedGSM_with_NonOverlapedInteraction = []
dataframe_GSM = pd.DataFrame()
dataframe_hmec =pd.DataFrame()

firstBin_have_promoter = False
secondBin_have_promoter = False


# def hash_start_end(num):
#     return int(int(num) // 5000) * 5000


try:
    # در این دیکشنری کلید را  چر , و ولیو را شروع  اچمک و شروع  اینترکشن قرار میدهیم که بعدا آنها را اسپلیت کنیم
    with open('h.csv', 'r') as hmecfile:
        csv_reader_hmec = csv.reader(hmecfile)
        dataframe_hmec = pd.DataFrame(csv_reader_hmec)

        for i in range(1, len(dataframe_hmec.index)):

            chr_of_hmec = int(dataframe_hmec.iloc[[i]][2].to_string(index=False))
            start_of_hmec = int(dataframe_hmec.iloc[[i]][3].to_string(index=False))
            end_of_hmec = int(dataframe_hmec.iloc[[i]][4].to_string(index=False))
            start_of_interaction_hmec = int(dataframe_hmec.iloc[[i]][7].to_string(index=False))
            end_of_interaction_hmec = int(dataframe_hmec.iloc[[i]][8].to_string(index=False))


            info_list = []

            for i in all_chars:
                if chr_of_hmec == i:
                    info_list.append(start_of_hmec)
                    info_list.append(end_of_hmec)
                    info_list.append(start_of_interaction_hmec)
                    info_list.append(end_of_interaction_hmec)
                    HMEC_dict[i].append(info_list)

                info_list = []

except IOError as e:
    print('Operation failed: %s' % e.strerror)

print(HMEC_dict)


with open('p.csv', 'r') as Promfile:

    csv_reader = csv.reader(Promfile)
    dataframe_promoter = pd.DataFrame(csv_reader)
    numberOf_allPromoters = (len(dataframe_promoter.index)) - 1

    for i in range(1, len(dataframe_promoter.index) ):
        begin_overflow = int(dataframe_promoter.iloc[[i]][1].to_string(index=False))
        end_overflow = int(dataframe_promoter.iloc[[i]][2].to_string(index=False))
        char_promoter = int(dataframe_promoter.iloc[[i]][5].to_string(index=False))

        if char_promoter in HMEC_dict.keys():  # char haye yeksan
            for i in range(len(HMEC_dict[char_promoter])):  # be tedad value haye key
                beggin_hmec = HMEC_dict[char_promoter][i][0]  # adad sotonaye start hmec
                end_hmec = HMEC_dict[char_promoter][i][1]
                beggin_interaction = HMEC_dict[char_promoter][i][2]
                end_interaction = HMEC_dict[char_promoter][i][3]

                if ( beggin_hmec <= begin_overflow ) and ( end_overflow <= end_hmec ) :
                    dataframe_hmec.loc[dataframe_hmec[3].eq(str(beggin_hmec)) & dataframe_hmec[7].eq(str(beggin_interaction)), 12] = True
                if ( beggin_interaction <= begin_overflow ) and ( end_overflow <= end_interaction ) :
                    dataframe_hmec.loc[dataframe_hmec[3].eq(str(beggin_hmec)) & dataframe_hmec[7].eq(str(beggin_interaction)), 13] = True


                if 10 <= ( beggin_hmec - end_overflow ) and  ( beggin_hmec - end_overflow )<=2010 :
                    dataframe_hmec.loc[dataframe_hmec[3].eq(str(beggin_hmec)) & dataframe_hmec[7].eq(str(beggin_interaction)), 12] = True
                if 10 <= ( beggin_interaction - end_overflow ) and  ( beggin_interaction - end_overflow )<=2010 :
                    dataframe_hmec.loc[dataframe_hmec[3].eq(str(beggin_hmec)) & dataframe_hmec[7].eq(str(beggin_interaction)), 13] = True


                if 10 <= (end_hmec - begin_overflow ) and (end_hmec - begin_overflow )<=2010 :
                    dataframe_hmec.loc[dataframe_hmec[3].eq(str(beggin_hmec)) & dataframe_hmec[7].eq(str(beggin_interaction)), 12] = True
                if 10 <= ( end_interaction - begin_overflow ) and ( end_interaction - begin_overflow )<=2010 :
                    dataframe_hmec.loc[dataframe_hmec[3].eq(str(beggin_hmec)) & dataframe_hmec[7].eq(str(beggin_interaction)), 13] = True





print(dataframe_hmec)

        # if start == end:
        #     if char_promoter in HMEC_dict.keys():  # char haye yeksan
        #         for i in range(len(HMEC_dict[char_promoter])): # be tedad value haye key
        #             beggin_hmec = HMEC_dict[char_promoter][i][0] # adad sotonaye start hmec
        #             beggin_interaction = HMEC_dict[char_promoter][i][1]
        #
        #             if start == beggin_hmec :
        #                 dataframe_hmec.loc[dataframe_hmec[3].eq( str(beggin_hmec) ) & dataframe_hmec[7].eq( str(beggin_interaction) ), 12] = True
        #
        #             if start == beggin_interaction:
        #                 dataframe_hmec.loc[dataframe_hmec[3].eq( str(beggin_hmec) ) & dataframe_hmec[7].eq( str(beggin_interaction) ), 13] = True
        #
        #
        # else:
        #     if char_promoter in HMEC_dict.keys():
        #         for i in range(len(HMEC_dict[char_promoter])): # be tedad value haye key
        #             beggin_hmec = HMEC_dict[char_promoter][i][0] #sotonaye hmec
        #             beggin_interaction = HMEC_dict[char_promoter][i][1]
        #             end_interaction = HMEC_dict[char_promoter][i][2]
        #
        #             if (beggin_hmec + 5000) - int(begin_overflow) >= 10:
        #                 dataframe_hmec.loc[dataframe_hmec[3].eq(str(beggin_hmec)) & dataframe_hmec[7].eq(str(beggin_interaction)), 12] = True
        #
        #             if (beggin_interaction + 5000) - int(begin_overflow) >= 10:
        #                 dataframe_hmec.loc[dataframe_hmec[3].eq(str(beggin_hmec)) & dataframe_hmec[7].eq(str(beggin_interaction)), 13] = True
        #
        #
        #             if int(begin_overflow) - beggin_hmec >=10 :
        #                 dataframe_hmec.loc[dataframe_hmec[3].eq(str(beggin_hmec)) & dataframe_hmec[7].eq(str(beggin_interaction)), 12] = True
        #
        #             if int(begin_overflow) - beggin_interaction >=10 :
        #                 dataframe_hmec.loc[dataframe_hmec[3].eq(str(beggin_hmec)) & dataframe_hmec[7].eq(str(beggin_interaction)), 13] = True
        #






        # if start == end:
        #     if char_promoter in HMEC_dict.keys():  # char haye yeksan
        #         for i in range(len(HMEC_dict[char_promoter])):
        #
        #             beggin_hmec = HMEC_dict[char_promoter][i][0]
        #             beggin_interaction = HMEC_dict[char_promoter][i][1]
        #
        #             if start == beggin_hmec and start == beggin_interaction:  # peyda kardane onaee k do sareshon poromotere
        #
        #                 if char_promoter not in deleted_dict.keys():
        #                     deleted_dict[char_promoter] = []
        #                     deleted_dict[char_promoter].append(HMEC_dict[char_promoter][i])
        #
        #                 else:
        #                     deleted_dict[char_promoter].append(HMEC_dict[char_promoter][i])
        #
        #
        # else:
        #     if char_promoter in HMEC_dict.keys():  # char haye yeksan
        #         for i in range(len(HMEC_dict[char_promoter])):
        #
        #             beggin_hmec = HMEC_dict[char_promoter][i][0]
        #             beggin_interaction = HMEC_dict[char_promoter][i][1]
        #
        #             if start == beggin_hmec and start == beggin_interaction:
        #                 if char_promoter not in deleted_dict.keys():
        #
        #                     deleted_dict[char_promoter] = []
        #                     deleted_dict[char_promoter].append(HMEC_dict[char_promoter][i])
        #
        #                 else:
        #                     deleted_dict[char_promoter].append(HMEC_dict[char_promoter][i])
        #
        #             if end == beggin_hmec and end == beggin_interaction:
        #                 if char_promoter not in deleted_dict.keys():
        #
        #                     deleted_dict[char_promoter] = []
        #                     deleted_dict[char_promoter].append(HMEC_dict[char_promoter][i])
        #
        #                 else:
        #                     deleted_dict[char_promoter].append(HMEC_dict[char_promoter][i])























































