'''

created on: 3rd July, 2015
author: Rohith Subrmanyam <rohithvsm@gmail.com>

'''

# required packages and modules
import os
import sys
import collections


def do(table):
    '''
    '''

    #sep = '^A'  # type this as control-v and control-a

    uniques_table = {} #collections.defaultdict(int)

    for tuple in table:
        key = tuple['id'] #'%s%s%s' % (tuple['type'], sep, tuple['id'])
        if key in uniques_table:
            if float(tuple['match_percentage']) > uniques_table[key][0]:
                uniques_table[key] = [float(tuple['match_percentage']), tuple['type'], tuple['alignment'] ]#max(uniques_table[key], tuple['match_percentage'])
        else:
            uniques_table[key] = [float(tuple['match_percentage']), tuple['type'], tuple['alignment'] ]

    result_table = []

    for key in uniques_table:
        uniq_tuple = {}
        #key_splits = key.split(sep)
        #type = key_splits[0]
        #id = key_splits[1]
        uniq_tuple['id'] = key #id
        uniq_tuple['type'] = uniques_table[key][1]#type
        uniq_tuple['match_percentage'] = uniques_table[key][0]
        uniq_tuple['alignment'] = uniques_table[key][2]
        result_table.append(uniq_tuple)

    return result_table


def main():
    args = sys.argv

    if len(args) > 1:
        alignment_filename = args[1]

        import load_alignments
        alignments = load_alignments.do(alignment_filename)

    else:
        alignments = [{'id':'P0', 'type':'r','match_percentage':34.88, 'alignment':'ABCABSGHYA'},
                {'id':'P0', 'type':'r','match_percentage':76.21, 'alignment':'HHGSQQTUASVLLMN'},
                {'id':'P1', 'type':'r','match_percentage':99.90, 'alignment':'JEDUIWDJ'},
                {'id':'P1', 'type':'r','match_percentage':99.90, 'alignment':'HGIWHEJHILWLD'}]

    print(do(alignments))


if __name__ == '__main__':
    main()
