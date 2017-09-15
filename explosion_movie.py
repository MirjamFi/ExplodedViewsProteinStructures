1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
31
32
33
34
35
36
37
38
39
40
41
42
43
44
45
46
47
48
49
50
51
52
53
54
55
56
57
58
59
60
61
62
63
64
65
66
67
68
69
70
71
72
73
74
75
76
77
78
79
80
81
82
83
84
85
86
87
88
89
90
91
92
93
94
95
96
97
98
99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
173
174
175
176
177
178
179
180
181
182
183
184
185
186
187
188
189
190
191
192
193
194
195
196
197
198
199
200
201
202
203
204
205
206
207
208
209
210
211
212
213
214
215
216
217
218
219
220
221
222
223
224
225
226
227
228
229
230
231
232
233
234
235
236
237
238
239
240
241
242
243
244
245
246
247
248
249
250
251
252
253
254
255
256
257
258
259
260
261
262
263
264
265
266
267
268
269
270
271
272
273
274
275
276
277
278
279
280
281
282
283
284
285
286
287
288
289
290
291
292
293
294
295
296
297
298
299
300
301
302
303
304
305
306
307
308
309
310
311
312
313
314
315
316
317
318
319
320
321
322
323
324
325
326
327
328
329
330
331
332
333
334
335
336
337
338
339
340
341
342
343
344
345
346
347
348
349
350
351
352
353
354
355
356
357
358
359
360
361
362
363
364
365
366
367
368
369
370
371
372
373
374
375
376
377
378
379
380
381
382
383
384
385
386
387
388
389
390
391
392
393
394
395
396
397
398
399
400
401
402
403
404
405
406
407
408
409
410
411
412
413
414
415
416
417
418
419
420
421
422
423
424
425
426
427
428
429
430
431
432
433
434
435
436
437
438
439
440
441
442
443
444
445
446
447
448
449
450
451
452
453
454
455
456
457
458
459
460
461
462
463
464
465
466
467
468
469
470
471
472
473
474
475
476
477
478
479
480
481
482
483
484
485
486
487
488
489
490
491
492
493
494
495
496
497
498
499
500
501
502
503
504
505
506
507
508
509
510
511
512
513
514
515
516
517
518
519
520
521
522
523
524
525
526
527
528
529
530
531
532
533
534
535
536
537
538
539
540
541
542
543
544
545
546
547
548
549
550
551
552
553
554
555
556
557
558
559
560
561
562
563
564
565
566
567
568
569
570
571
572
573
574
575
576
577
578
579
580
581
582
583
584
585
586
587
588
589
590
591
592
593
594
595
596
597
598
599
600
601
602
603
604
605
606
607
608
609
610
611
612
613
614
615
616
617
618
619
620
621
622
623
624
625
626
627
628
629
630
631
632
633
634
635
636
637
638
639
640
641
642
643
644
645
646
647
648
649
650
651
652
653
654
655
656
657
658
659
660
661
662
663
664
665
666
667
668
669
670
671
672
673
674
675
676
677
678
679
680
681
682
683
684
685
686
687
688
689
690
691
692
693
694
695
696
697
698
699
700
701
702
703
704
705
706
707
708
709
710
711
712
713
714
715
716
717
718
719
720
721
722
723
724
725
726
727
728
729
730
731
732
733
734
735
736
737
738
739
740
741
742
743
744
745
746
747
748
749
750
751
752
753
754
755
756
757
758
759
760
761
762
763
764
765
766
767
768
769
770
771
772
773
774
775
776
777
778
779
780
781
782
783
784
785
786
787
788
789
790
791
792
793
794
795
796
797
798
799
800
801
802
803
804
805
806
807
808
809
810
811
812
813
814
815
816
817
818
819
820
821
822
823
824
825
826
827
828
829
830
831
832
833
834
835
836
837
838
839
840
841
842
843
844
845
846
847
848
849
850
851
852
853
854
855
856
857
858
859
860
861
862
863
864
865
866
867
868
869
870
871
872
873
874
875
876
877
878
879
880
881
882
883
884
885
886
887
888
889
890
891
892
893
894
895
896
897
898
899
900
901
902
903
904
905
906
907
908
909
910
911
912
913
914
915
916
917
918
919
920
921
922
923
924
925
926
927
928
929
930
931
932
933
934
935
936
937
938
939
940
941
942
943
944
945
946
947
948
949
950
951
952
953
954
955
956
957
958
959
960
961
962
963
964
965
966
967
968
969
970
971
972
973
974
975
976
977
978
979
980
981
982
983
984
985
986
987
988
989
990
991
992
993
994
995
996
997
998
999
1000
1001
1002
1003
1004
1005
1006
1007
1008
1009
1010
1011
1012
1013
1014
1015
1016
1017
1018
1019
1020
1021
1022
1023
1024
1025
1026
1027
1028
1029
1030
1031
1032
1033
1034
1035
1036
1037
1038
1039
1040
1041
1042
1043
1044
1045
1046
1047
1048
1049
1050
1051
1052
1053
1054
1055
1056
1057
1058
1059
1060
1061
1062
1063
1064
1065
1066
1067
1068
1069
1070
1071
1072
1073
1074
1075
1076
1077
1078
1079
1080
1081
1082
1083
1084
1085
1086
1087
1088
1089
1090
1091
1092
1093
1094
1095
1096
1097
1098
1099
1100
1101
1102
1103
1104
1105
1106
1107
1108
1109
1110
1111
1112
1113
1114
1115
1116
1117
1118
1119
1120
1121
1122
1123
1124
1125
1126
1127
1128
1129
1130
1131
1132
1133
1134
1135
1136
1137
1138
1139
1140
1141
1142
1143
1144
1145
1146
1147
1148
1149
1150
1151
1152
1153
1154
1155
1156
1157
1158
1159
1160
1161
1162
1163
1164
1165
1166
1167
1168
1169
1170
1171
1172
1173
1174
1175
1176
1177
1178
1179
1180
1181
1182
1183
1184
1185
1186
1187
1188
1189
1190
1191
1192
1193
1194
1195
1196
1197
1198
1199
1200
1201
1202
1203
1204
1205
1206
1207
1208
1209
1210
1211
1212
1213
1214
1215
1216
1217
1218
1219
1220
1221
1222
1223
1224
1225
1226
1227
1228
1229
1230
1231
1232
1233
1234
1235
1236
1237
1238
1239
1240
1241
1242
1243
1244
1245
1246
1247
1248
1249
1250
1251
1252
1253
1254
1255
1256
1257
1258
1259
1260
1261
1262
1263
1264
1265
1266
1267
1268
1269
1270
1271
'''
DESCRIPTION
 
explosion creates a movie of the exploded view of a molecule. 
If there are multiple objects given for explosion, first they get spatially 
separated and then explode individually one after the other. The order of 
explosion is the order of given list of objects.
There are two types of explosion direction:
- 'com' (default):  the centers of mass (com) of the chains of the object to be 
                    tranlated and the object are calulated and the single chains 
                    are translated along a vector through the chain's com and 
                    object's com.
- 'canonical':      the dimensions of a box around the object are used to select 
                    the two longest edges and so the axes to translate along in 
                    a consistent distance
                     
If only a part of the object shall be translated the object can be given as 
complex to make sure the part is not translated into the object.
 
AUTHOR
 
    Mirjam Figaschewski
    mirjam_figaschewski (at) web.de
 
DEPENDENCIES:
get_colors.py (https://pymolwiki.org/index.php/Get_colors) and 
center_of_mass.py (https://pymolwiki.org/index.php/Center_of_mass)
viewpoints.py (https://github.com/julianheinrich/viewpoints)
 
EXAMPLES
 
explosion 5gmz
'''
 
from pymol import cmd 
import math
import center_of_mass as cenma ## calulate center of mass
from pymol import stored ## import stored for passing data back and forth
import get_colors ## get random colors
import time
from operator import itemgetter
from viewpoints import best_view
import os.path
 
def initialize_movie(selected = None, frames = "100"):
    ''' DESCRIPTION: 
        initial setup of a movie'''
    cmd.set('matrix_mode', 1)
    cmd.set('movie_panel', 1)
    cmd.set('scene_buttons', 1)
    cmd.set('cache_frames', 1)
    cmd.config_mouse('three_button_motions', 1)
    # cmd.set('movie_panel', 0) ## hide movie panel
    cmd.set('movie_panel_row_height', 1)
    cmd.set('movie_fps', 10)
     
    cmd.mset('1 x' + frames)
    if selected:
        cmd.orient(selected)
 
def remove_solvents(exclude = "", cutoff =10, storedLigands = None):
    cutoff = int(cutoff)
    '''DESCRIPTION:
        removes ligands/solvents with occurence higher than cutoff in cc-counts.tdd
        ligands in excluded are not removed
    '''
     
    if not os.path.exists('./cc-counts.tdd'):
        import urllib
        urllib.urlretrieve ("http://ligand-expo.rcsb.org/dictionaries/cc-counts.tdd", 
                                "cc_counts.tdd")
    if not os.path.exists('./cc-counts.tdd'):
        sys.exit('Could not download http://ligand-expo.rcsb.org/dictionaries/cc-counts.tdd, please try manually.')
     
    solvents = {}
    aa = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
    cmd.remove('solvent')
     
    import csv
    removedLigands = []
    with open('./cc-counts.tdd') as csvfile:
        countreader = csv.reader(csvfile, delimiter='\t')
        for row in countreader:
            if row[0] != 'id':
                ligand = row[0]
                ligcount = int(row[1])
                if ligcount > cutoff and not ligand in aa:
                    if not exclude or ligand not in exclude:
                        cmd.remove('resn ' + ligand)
                        if storedLigands:
                            for lig in storedLigands:
                                if lig.split('-')[0] == ligand:
                                    removedLigands.append(lig)
 
    for removelig in removedLigands:
        storedLigands.discard(removelig)
 
def get_ligands():
    ''' DESCRIPTION:
        get all ligands (organic) from complex, return them as a set'''
    stored.ligands = []
 
    ## create the original selection of all organic atoms 
    cmd.select('allOrg', 'org')
 
    ##' temporary selection
    cmd.select('temp', 'none')
 
    ## while initial selection is not empty, 'pop' from it and query the atom's 
    ## molecule's name
    while cmd.count_atoms('allOrg') != 0:
 
        ## pop--peek, rather
        cmd.select('temp', 'first allOrg')
         
        ## store the residue name
        cmd.iterate('temp', 'stored.ligands.append("%s-%s" %(resn,resi))')
         
        ## remove by molecule
        cmd.select('allOrg', 'allOrg and not (bm. temp)')
         
    cmd.delete('allOrg')
    cmd.delete('temp')
 
    return set(stored.ligands)
     
def calc_COM(transname):
    '''DESCRIPTION:
        calculate CenterOfMass of an object 
    '''
    cenma.com(transname, state = 1)
    cmd.zoom(transname + '_COM')
    pos = cmd.get_position(transname + '_COM')
    cmd.delete(transname + '_COM')
    cmd.zoom('all', complete = 1)
    cmd.mview('store')
 
    return pos
 
def calcTransFac(sele):
    '''DESCRIPTION 
        calculate translation factor according to size of selection
    '''
    protDim_min, protDim_max = cmd.get_extent(sele)
    transFac = math.sqrt((protDim_max[0] - protDim_min[0]) ** 2 +
                     (protDim_max[1] - protDim_min[1]) ** 2 +
                     (protDim_max[2] - protDim_min[2]) ** 2)/4
    return transFac
 
def transAxes(selected):
    '''DESCRIPTION:
        get dimensions of structure to decide along which axes to translate
    '''
    box1_min, box1_max = cmd.get_extent(selected)
 
    x_axis = abs(box1_min[0] - box1_max[0])
    y_axis = abs(box1_min[1] - box1_max[1])
    z_axis = abs(box1_min[2] - box1_max[2])
 
    axes = {"x_axis" : x_axis, "y_axis":y_axis, "z_axis":z_axis}
    axes = sorted(axes.items(), key=itemgetter(1), reverse = True)
 
    transFac = calcTransFac(selected)
    transVec = [0,0,0]
    if axes[0][0] == "x_axis" or axes[1][0] == "x_axis":
        transVec[0] = transFac
    if axes[0][0] == "y_axis" or axes[1][0] == "y_axis":
        transVec[1] = transFac  
    if axes[0][0] == "z_axis" or axes[1][0] == "z_axis":
        transVec[2] = transFac
         
    return transVec
     
def create_objects(chains, selected, storedLigands, chainsCOMS, complexXYZ, dim, 
                    colorBinding, chainAndLigand = None, typeOfExplosion = 'com'):
    ''' DESCRIPTION:
        for every chain and ligand a unique object is created, coloring is 
        applied and labels created
    '''
 
    label_objec = {}
     
    ## names of chains
    cNames = []
     
    ## chains and their labels
    chainAndLabel = {}
     
    ## chains and according ligands
    if not chainAndLigand:
        chainAndLigand = {}
     
    ## ligands and according chains
    ligandAndChain={}
     
    ## COMS of ligands/chains
    if typeOfExplosion == 'com':
        ligandsCOMS = {}
        chainsCOMS = {}
         
    ## store coloring
    objColor = {}
     
    for c in chains:    
        if colorBinding == 'chain':
            ## color each chain individually (red is color for binding site)
            col = get_colors.get_random_color()
            while col in objColor.values() or \
                    col in ('red','white', 'black', 'marine'):
                col = get_colors.get_random_color()
            cmd.color(col, selected + '& chain ' + c) 
             
        if colorBinding == 'contact':
            ## color all chains in gray
            cmd.color('gray', selected + '& chain ' + c)
                     
        ## create object for chain with name
        chainname = selected + '_chain' + c
        cNames.append(chainname)
        cmd.extract(chainname, selected + ' & chain '+c)
        if colorBinding == 'chain':
            objColor.update({chainname:col})
 
        if typeOfExplosion=='com':
             
            ## calculate and save COM of chain
            chainsCOMS[chainname] = calc_COM(chainname)
            chainCOM = chainsCOMS[chainname]
             
            ## label chains
            labels = label_obj(chainname, chainCOM, complexXYZ, dim, chainAndLabel)
            if len(labels) > 1:
                chainAndLabel = labels[0]
                label_objects_new = labels[1]
            elif len(labels) == 1:
                label_objects_new = labels
             
            if label_objects_new:
                label_objec.update(label_objects_new)
 
        ## store object for movie
        f = 1
        if typeOfExplosion == 'com':
            cmd.frame(f)
            cmd.mview('store', object=chainname)
            cmd.mview('store', object="label" + chainname)
            cmd.zoom('all', complete = 1)
            cmd.show('dashes')
         
        ## create objects for ligands
        if storedLigands:
            i = 1
            for l in storedLigands:
                resname, resid = l.split('-')
             
                ## select ligand on chain   
                selection = chainname + ' & resn ' + resname + "& resi " + resid
 
                if cmd.count_atoms(selection) > 0:
                    if typeOfExplosion == 'com':
                        label_objec_new = get_ligand_chain_pair(l,i, selection, 
                                                    chainname, chainAndLigand, 
                                                    ligandAndChain, complexXYZ, dim, 
                                                    label_objec, 'com',ligandsCOMS)
                        if label_objects_new:
                            label_objec.update(label_objects_new)
                    else:
                        get_ligand_chain_pair(l, i,selection, chainname, 
                                            chainAndLigand, ligandAndChain, 
                                            complexXYZ, dim, label_objec, 
                                            typeOfExplosion = 'canonical')
                         
                    cmd.frame(f)
                    cmd.mview('store', object = 'all')
                    i += 1
                else:
                    continue
            ## if there is no ligand on chain, just keep chain and its label
            if typeOfExplosion == 'com':
                noLigands(chainname, chainAndLigand, label_objec)
            else:
                noLigands(chainname, chainAndLigand)
        ## no ligands in structure
        else:
            if typeOfExplosion == 'com':
                noLigands(chainname, chainAndLigand, label_objec)
            else:
                noLigands(chainname, chainAndLigand)
                 
    if typeOfExplosion == 'com':    
        return cNames, chainAndLabel, chainAndLigand, ligandAndChain, ligandsCOMS,chainsCOMS, f, label_objec, objColor
    else:
        return cNames, chainAndLabel, chainAndLigand, ligandAndChain, f, label_objec, objColor
 
def get_ligand_chain_pair(l, i, selection, chainname, chainAndLigand, 
                            ligandAndChain, complexXYZ, dim, label_objects, 
                            typeOfExplosion = 'com', ligandsCOMS = None):
    '''DESCRIPTION:
        get pair of ligand and respective chain'''
 
    ## create object for ligand
     
    ligandname = chainname + '_' + l.split('-')[0] +'_'+ str(i)
    cmd.extract(ligandname, selection)
 
    ligandAndChain[ligandname] = chainname
     
    if typeOfExplosion == 'com':
        ## calculate COM of ligand
        ligandsCOMS[ligandname] = calc_COM(ligandname)
        ligandCOM = ligandsCOMS[ligandname]
        ## label ligands            
        label_objects_new = label_obj(ligandname, ligandCOM, complexXYZ, dim)
     
        ## group chains and respective ligands to translate  
        ## them together at first (including their labels)
        cmd.group(ligandname + "_", 
                        chainname + " " + ligandname + " " + "label" + \
                        chainname + " " + "label" + ligandname + " " + \
                        label_objects_new[ligandname][0] + " " + \
                        label_objects[chainname][0])
        chainAndLigand[chainname].append(ligandname + "_")
 
        label_objects.update(label_objects_new)
         
        return label_objects
     
    if typeOfExplosion == 'canonical':
        ## group chains and respective ligands to translate  
        ## them together at first (including their labels)
        cmd.group(ligandname + "_", chainname + " " + ligandname)
        chainAndLigand[chainname].append(chainname + "_" + l.split('-')[0] +\
                                                '_'+ str(i) + "_")          
         
def calc_label_positions_circular(ch,chainCOM, complexXYZ, dim):
    ''' DESCRIPTION:
        set labels on sphere around complex object with connection line to com
        of according chain/ligand
    '''
     
    ## sphere
    min = dim[0]
    max = dim[1]
    R = math.sqrt((min[0]-max[0])**2 + (min[1]-max[1])**2 + (min[2]-max[2])**2) #euclidean
    R = R/2
     
    P1 = chainCOM
    P2 = complexXYZ
    x0 = P1[0]
    y0 = P1[1]
    z0 = P1[2]
    x1 = P2[0]
    y1 = P2[1]
    z1 = P2[2]
    dx = x1 - x0
    dy = y1 - y0
    dz = z1 -  z0 
    cx = P2[0]
    cy = P2[1]
    cz = P2[2]
     
    ## intersection point sphere - com-vector
    a = dx*dx + dy*dy + dz*dz
    b = 2*dx*(x0-cx)+2*dy*(y0-cy)+2*dz*(z0-cz)
    c = cx*cx+cy*cy+cz*cz+x0*x0+y0*y0+z0*z0+-2*(cx*x0+cy*y0+cz*z0)-R*R
     
    discriminant = b**2 - 4*a*c
     
    t = (-b-math.sqrt(discriminant))/(2*a)
     
    x = x0 + t*dx                 
    y = y0 + t*dy                 
    z = z0 + t*dz
     
    n1 = ''.join(map(str,P1))
 
    cmd.pseudoatom('_pt1' + n1, pos=P1)
    cmd.pseudoatom('_pt2' + n1, pos=[x,y,z])
    cmd.hide('nonbonded', "_pt1" + n1)
    cmd.hide('nonbonded', "_pt2" + n1)
    cmd.distance('_'+ch+'_label', '_pt1' + n1, '_pt2' + n1)
    f = 1
    cmd.frame(f)
    cmd.hide('labels', '_'+ch+'_label')
    cmd.mview('store', object='_'+ch+'_label')
    return [x,y,z], ['_'+ch+'_label']
         
def calc_label_position_flush(chains, transVec, f, chainAndLabel=None):
    ''' DESCRIPTION:
        set labels on left or right (up/down) of object, connection line to com
        of according chain/ligand
    '''
    chainsComs = {}
    label_objects = {}
 
    ## Determine the positions of anchor points,
    size = {}
    for ch in chains:
        chainsComs[ch] = calc_COM(ch)
        size[ch] = cmd.count_atoms(ch)
    maxChain = max(size, key=size.get)
     
    maxdim = cmd.get_extent(maxChain)
    cmd.pseudoatom('_min'+maxChain, pos=maxdim[0])
    cmd.hide('nonbonded', '_min'+maxChain)
    cmd.pseudoatom('_max'+maxChain, pos=maxdim[1])
    cmd.hide('nonbonded', '_max'+maxChain)
    labellength = cmd.get_distance('_min'+maxChain,'_max'+maxChain)
     
    ## Choose a pivot point
    mean_xpos = sum(chainsComs[ch][0]for ch in chainsComs.keys())/len(chainsComs)
     
    mean_ypos = sum(chainsComs[ch][1] for ch in chainsComs.keys())/len(chainsComs)
     
    mean_zpos = sum(chainsComs[ch][2] for ch in chainsComs.keys())/len(chainsComs)
     
    cmd.pseudoatom('_mean', pos=[mean_xpos,mean_ypos,mean_zpos])
    cmd.hide('nonbonded', '_mean')
 
    ## Assign labels to the left and right region by comparing
    ## their anchors with the pivot point.
    for ch in chainsComs:
        cmd.pseudoatom('_com_' + ch, pos = chainsComs[ch])
        cmd.hide('nonbonded', '_com_' + ch)
         
        if transVec[0] > 0 and transVec[1]:
            x_pos = chainsComs[ch][0]
            y_pos = chainsComs[ch][1]
            if chainsComs[ch][2] > mean_zpos:
                z_pos = chainsComs[ch][2] + labellength
            else:
                z_pos = chainsComs[ch][2] - labellength
             
        elif transVec[1] > 0 and transVec[2] > 0:
            y_pos = chainsComs[ch][1]
            z_pos = chainsComs[ch][2]
            if chainsComs[ch][0] > mean_xpos:
                x_pos = chainsComs[ch][0] + labellength
            else:
                x_pos = chainsComs[ch][0] - labellength
                 
        elif transVec[2] > 0 and transVec[0] > 0:
            x_pos = chainsComs[ch][0]
            z_pos = chainsComs[ch][2]
            if chainsComs[ch][1] > mean_ypos:
                y_pos = chainsComs[ch][1] + labellength
            else:
                y_pos = chainsComs[ch][1] - labellength
                 
        cmd.pseudoatom('_'+ch+'_pos', pos = [x_pos, y_pos, z_pos])
        cmd.distance('_'+ch + '_label','_'+ch+'_pos', '_com_' + ch)
        cmd.hide('nonbonded', '_'+ch+'_pos')    
        cmd.frame(f)
        cmd.pseudoatom('label' + ch, pos = [x_pos, y_pos, z_pos])
        if len(ch.split('_')) > 2:
            cmd.label("label" + ch, "'%s'" %'_'.join(ch.split('_')[-2:]))
        else:
            cmd.label("label" + ch, "'%s'" %ch.split('_')[-1])
        cmd.hide('labels', '_'+ch + '_label')
        cmd.hide('nonbonded', 'label' + ch)
         
        chainAndLabel[ch] = "label" + ch
        label_objects.update({ch : '_'+ch + '_label'})
    return [chainAndLabel, label_objects]
 
def label_obj(chainname, chainCOM, complexXYZ, dim, chainAndLabel = None):
    '''DESCRIPTION:
        create an label for given chain (if chainAndLabel set) or ligand  
    '''
    ## create pseudoatom for object to be labeled
    label_pos, lab_obj = calc_label_positions_circular(chainname, chainCOM, 
                                                                complexXYZ, dim)
     
    cmd.pseudoatom("label" + chainname, pos=label_pos)
     
    ## label object
    if not chainAndLabel:
        if len(chainname.split('_')) > 2:
            ligname = '_'.join(chainname.split("_")[-2:])
            cmd.label("label" + chainname, "'%s'" %ligname)
        else:
            cmd.label("label" + chainname, "'%s'" %chainname.split('_')[-1])
    else:
        cmd.label("label" + chainname, "'%s'" %chainname)
     
    ## save label of chain in dictonary
    if chainAndLabel:
        chainAndLabel[chainname] = "label" + chainname
 
    ## store ligand label in movie
    if not chainAndLabel:
        cmd.mview('store', object = "label" + chainname)
         
    ## hide pseudoatom represantation
    cmd.hide('nonbonded', "label" + chainname)
     
    ## return dictionary of chain and its label
    if chainAndLabel:
        return [chainAndLabel, {chainname : lab_obj}]
    else:
        return {chainname : lab_obj}
 
def noLigands(chainname, chainAndLigand, label_objec = None):
    ''' DESCRIPTION:
        create group for chain which do not have a ligand
    '''
    if chainname in chainAndLigand.keys() and not chainAndLigand[chainname]:
        del chainAndLigand[chainname]
    if label_objec:
        cmd.group(chainname + "_", chainname + " " + "label" + chainname + " "\
                    + label_objec[chainname][0])
    else:
        cmd.group(chainname + "_", chainname + " " + "label" + chainname)
         
def color_binding(chainname, ligandname, contactColors = None):
    ''' DESCRIPTION:
        color binding site in red'''
    binding = 'byres (' + chainname + ' nto. 3.6 of ' + ligandname +')'
    cmd.select('_inter1', binding)
     
    binding2 = 'byres (' + ligandname + ' nto. 3.6 of ' + chainname +')'
    cmd.select('_inter2', binding2)
     
    if not contactColors:
        cmd.color('red', '_inter1')
         
    if contactColors:
        col = get_colors.get_random_color()
        while col in contactColors or col in ('white', 'black', 'marine'):
            col = get_colors.get_random_color()
        cmd.color(col, '_inter1') 
        cmd.color(col, '_inter2')
        contactColors.append(col)   
    cmd.delete('_inter')
 
def color_contact(cNames, objColor, colorBinding):
    ''' DESCRIPTION:
        color contact site between two chains by color of the other chain
    '''
    if colorBinding == 'contact':
        contactColors = ['gray']
    i = 1
    for chain in cNames:
        for ch in cNames[i:]:
            if ch != chain:
             
                binding = 'byres (' + chain + ' nto. 5 of ' + ch +')'
                cmd.select('_contact1', binding)
                 
                binding = 'byres (' + ch + ' nto. 5 of ' + chain +')'
                cmd.select('_contact2', binding)
                 
                if colorBinding == 'chain':
                    cmd.color(objColor[ch], '_contact1')
                    cmd.color(objColor[chain], '_contact2')
                     
                if colorBinding == 'contact':
                    col = get_colors.get_random_color()
                    while col in contactColors or \
                                    col in ('white', 'black', 'marine'):
                        col = get_colors.get_random_color()
                    cmd.color(col, '_contact1') 
                    cmd.color(col, '_contact2') 
                    contactColors.append(col)           
        i+=1
    if colorBinding == 'contact':
        return contactColors
 
def create_complex(chains, obj):
    '''DESCRIPTION:
        for list of chains create object from chains belonging to obj
    '''
    sel = ''
    for ch in chains[obj]:
        sel = sel + obj+'_chain' + ch + ' '
    cmd.select('_sel' + obj, sel)
     
def isColliding(sel1, sel2):
    ''' DESCRIPTION:
        calculate if two bounding boxes are colliding (boxes are axis-aligned)
    '''
    box1_min, box1_max = cmd.get_extent(sel1)
    box2_min, box2_max = cmd.get_extent(sel2)
 
    ## min, max vertices of box 1
    x_min1 = box1_min[0]
    y_min1 = box1_min[1]
    z_min1 = box1_min[2]
    x_max1 = box1_max[0]
    y_max1 = box1_max[1]
    z_max1 = box1_max[2]
     
    ## min, max vertices of box 2
    x_min2 = box2_min[0]
    y_min2 = box2_min[1]
    z_min2 = box2_min[2]
    x_max2 = box2_max[0]
    y_max2 = box2_max[1]
    z_max2 = box2_max[2]
     
    ## test if any of the vertices of box 1 are in any of the faces of box 2 
    ## and vice versa
    isColliding =   (x_min1 <= x_max2 and x_max1 >= x_min2) and \
                    (y_min1 <= y_max2 and y_max1 >= y_min2) and \
                    (z_min1 <= z_max2 and z_max1 >= z_min2)
 
    return isColliding
 
def best_view_objects(label = False):
    '''DESCRIPTION: 
        return string of all objects condsidered for best_view calulation 
    '''
    view_objects = " "
    for obj in cmd.get_names('objects'):
            if not label:
                if not obj.startswith('_') and not obj.startswith('label'):
                    view_objects += obj + " "
            else:
                if not obj.startswith('_'):
                    view_objects += obj + " "
    return view_objects
                 
def translate_selection(originXYZ, transXYZ, transname, factor = 1, f = 1, 
                                                                group = None):
    ''' DESCRIPTION:
        translate an object relative to complex of origin using center of mass
    '''
    ## vector between COMs to translate chain
    dist=math.sqrt((transXYZ[0]-originXYZ[0])**2 +
                    (transXYZ[1]-originXYZ[1])**2 +
                    (transXYZ[2]-originXYZ[2])**2)
 
    vector=((transXYZ[0]-originXYZ[0])/dist, 
            (transXYZ[1]-originXYZ[1])/dist, 
            (transXYZ[2]-originXYZ[2])/dist)
    trans_vec = [x * factor for x in vector] 
 
    ## translate chain with vector
    cmd.frame(f)
    cmd.translate(trans_vec, object=transname, camera=0)
    if group:
        cmd.mview('store', object=group)
        cmd.mview('interpolate', object=group)
    else:
        cmd.mview('store', object=transname)
        cmd.mview('interpolate', object=transname)
     
    store_view(group = True, all = True)
     
def com_translation(cname, chainAndLigand, complexXYZ, chainXYZ, transFac, f):
    '''DESCRIPTION:
        translate chain via COM
    '''    
    translate_selection(complexXYZ, chainXYZ, cname + '_', transFac, f, cname + "_")
    ## if chain contains ligand, translate ligand also
    if cname in chainAndLigand.keys():
         
        for l in chainAndLigand[cname]:
            translate_selection(complexXYZ, chainXYZ, l, transFac, f, l)
 
def canonical_translation(ch, i, transVec, chainAndLigand, label_objects):
    '''DESCRIPTION:
        translate chain ch canonical
    '''
    cmd.translate([x * i for x in transVec], object=ch + "_", camera=0)
    cmd.mview('interpolate', object = ch+"_")
    cmd.mview('store', object = ch+"_")
     
    ## if chain contains ligand, translate ligand also
    if ch in chainAndLigand.keys():
        for l in chainAndLigand[ch]:
            cmd.translate([x * i for x in transVec], object=l, camera=0)
            cmd.mview('interpolate', object=l)
            cmd.mview('store', object = l)  
         
def store_view(obj = None, group = False, all = True):
    ''' DESCRIPTION:
        store selected view 
    '''
         
    ## store all objects
    if group:
        for group in cmd.get_names("objects"):
                cmd.mview('store', object = group)
     
    ## zoom out on everything
    if all:
        cmd.zoom('all', complete=1)
     
    ## store single object
    if obj:
        cmd.mview('store', object = obj)
         
    ## store camera position
    cmd.mview('store')
 
def show_labels(label_objects):
    cmd.show('labels')
    for l in label_objects.values():
        if type(l) is list:
            cmd.hide('labels', l[0])
        else:
            cmd.hide('labels', l)
    cmd.show('dashes')
     
def com_explosion(selected, label_objects, cNames, chainAndLigand, 
                    ligandAndChain, chainsCOMS, ligandsCOMS, transFac, 
                    complex = None, com = None, frame = 1):
    ''' DESCRIPTION:
        create a movie for an exploded view of given protein. Explosion along
        vector of complex' COM and selection's COM.
    '''
    start_time = time.clock()
     
    f = frame 
             
    ''' Explosion: Create an object of every chain and for related ligands. 
                Calulate COM for complex and according chains and 
                translate chain along vector between COMs of complex and chain. 
                For every chain do the same with its ligands.
    '''
         
    ## COM of complex
    ## if only some chains are selected and source complex is given,
    ## select source complex for COM calculation
    if com:
        complexXYZ = com
    if not com:
        complexXYZ = calc_COM(selected)
    f = f + 10
    ''' translate chains '''
    for chainname in cNames:
        ## COM of chain
        chainXYZ = chainsCOMS[chainname]
        ## only translate chain if there are multiple chains
        if len(cNames) > 1:  
            ## translate chains and ligands
            for cname in cNames:
                ## check if current chain is colliding with any other chain
                if cname != chainname:
                    while isColliding(chainname, cname):
                        ## if collision, translate all chains in selection
                        for cname in cNames:
                            cXYZ = chainsCOMS[cname]
                            com_translation(cname, chainAndLigand, complexXYZ,
                                                cXYZ, transFac, f)
                         
            store_view(group = True, all = True)
                                     
        ## if only one chain is selected, translate it and its ligand
        else:
            cmd.frame(f)
            cmd.zoom(selected)
            store_view(chainname + '_')
             
            ## if source complex is known use its COM to translate
            if complex:
                while isColliding(chainname, complex):
                    com_translation(chainname, chainAndLigand, complexXYZ,
                                                chainXYZ, transFac, f)
            else:
                if chainname in chainAndLigand.keys():
                        cmd.translate([transFac, transFac, transFac], 
                                        object = chainAndLigand[chainname])
                else:
                    cmd.translate([transFac, transFac, transFac], 
                                        object = chainname + "_")
                                         
            store_view(chainname + '_')
 
    ## store objects for movie
    f = f + 20
    cmd.frame(f)
    if len(cNames) == 1:
        cmd.zoom('all', complete=1)
    store_view(group=True, all = True)
 
    ''' translate ligands '''
    if ligandAndChain:
        f = f + 15
        cmd.frame(f)
        for ligand in ligandAndChain.keys():
            cXYZ = chainsCOMS[ligandAndChain[ligand]]
            ligandXYZ = ligandsCOMS[ligand]
             
            if complex:
                condition = isColliding(ligand, complex) \
                            or isColliding(ligand, ligandAndChain[ligand])
            else:
                condition = isColliding(ligand, ligandAndChain[ligand])
            while condition:
                translate_selection(cXYZ, ligandXYZ, ligand, transFac/2, f, 
                                    chainAndLigand[ligandAndChain[ligand]])
                translate_selection(cXYZ, ligandXYZ, "label" + ligand, transFac/2, 
                                    f, chainAndLigand[ligandAndChain[ligand]])
                translate_selection(cXYZ, ligandXYZ, label_objects[ligand], 
                                    transFac/2, f)
                translate_selection(cXYZ, ligandXYZ, label_objects[ligand][0], 
                                    transFac/2, f)
                if complex:
                    condition = isColliding(ligand, complex) \
                            or isColliding(ligand, ligandAndChain[ligand])
                else:
                    condition = isColliding(ligand, ligandAndChain[ligand])
         
        f = f + 30
        ## show labels
        cmd.frame(f)
        view_objects = " "
        store_view(group=True, all = True)
    return f
 
def canonical_explosion(selected, label_objects, cNames, chainAndLigand, transVec, 
                            ligandAndChain = None, complex = None, frame = 1):
    '''DESCRIPTION: translate a selection in canonical direction and create
        a movie from it'''
    start_time = time.clock()
    cmd.orient(selected)
     
    f = frame
    ''' translate chains'''
    i = 1
    cmd.frame(f)
    if complex:
        for chain in cNames:
            while isColliding(chain, complex):
                for ch in cNames:
                    canonical_translation(ch, i, transVec, chainAndLigand, \
                                            label_objects)
                    i += 1
    i = 1
    for chain in cNames:
        for c in cNames:
            if c != chain:
                while isColliding(chain, c):
                    for ch in cNames[0:]:
                        canonical_translation(ch, i, transVec, chainAndLigand,\
                                                label_objects)
                        i += 1
     
    store_view(group = True, all = True)
     
    f = f + 30
    cmd.frame(f)
    store_view(group=True, all = True)
         
    ''' translate ligands'''
    if ligandAndChain:
        f = f + 15
        cmd.frame(f)
        for ligand in ligandAndChain.keys():
            if complex:
                condition = isColliding(ligand, complex) \
                            or isColliding(ligand, ligandAndChain[ligand])
            else:
                condition = isColliding(ligand, ligandAndChain[ligand])
            while condition:
                cmd.translate([x * 1/4  for x in transVec], object=ligand)
                cmd.translate([x * 1/4 for x in transVec], 
                                                    object="label" + ligand)
                cmd.translate([x * 1/4 for x in transVec], 
                                                    object='_'+ligand+'_label')
                if complex:
                    condition = isColliding(ligand, complex) \
                            or isColliding(ligand, ligandAndChain[ligand])
                else:
                    condition = isColliding(ligand, ligandAndChain[ligand])
            store_view(group = True, all = True)
         
        f = f + 30
        cmd.frame(f)
        store_view(group=True, all = True)
 
    return f
     
def explosion(selected = ' ', typeOfExplosion = 'com', complex = None, 
                removeSolvents = True, exclude = None, cutoff = 10, 
                colorBinding = 'contact', showlabels = False):
    '''DESCRIPTION:
        perform an explosion of selected object(s) given in a list and create   
        a movie. If two objects are given they will be separateted and then 
        exploded individually.
    '''
    selected = selected.split()
    if exclude:
        exclude = exclude.split()
    if not typeOfExplosion in ['com', 'canonical']:
        sys.exit("Specify explosion: com (default) or canonical.") 
         
    if not selected:
        sys.exit("Please give a list of selected objects (or just one) to translate.") 
         
    if complex:
        if cmd.count_atoms(complex) == 0:
            complex = None
            print 'No need of complex'
     
    if not colorBinding in ['chain', 'contact', 'none']:
        sys.exit("Specify a color for chains: none (do not color), chain, or contact.")
             
    '''setup of selected sturcture'''
    ## case sensitive for chain ids
    cmd.set('ignore_case', 'off') 
    # cmd.set('ray_trace_mode', 1)
    # cmd.set('ray_trace_frames', 'on')
     
    cmd.set('dash_color', 'marine')
    cmd.set('dash_round_ends', 'off')
    cmd.set('dash_width', 1)
    cmd.set('dash_gap', 0)
         
    '''preparation and translation'''
    if len(selected) >= 1:
        ## initialize movie
        start_time = time.clock()
        if len(selected) > 1:
            initialize_movie(frames = str(100*len(selected)))
        else:
            initialize_movie(frames = str(130))
         
        ## get ligands
        storedLigands = get_ligands()
         
        ## remove solvents
        if removeSolvents == True:
            if storedLigands:
                remove_solvents(exclude, cutoff, storedLigands)
            else:
                remove_solvents(exclude, cutoff)
         
        ## lists for further computations
        s = ''
        chains = {}
        cNames = []
        chainAndLabel = {}
        chainAndLigand = {}
        ligandAndChain = {}
        chainsAndObj = {}
        label_objects = {}
        objColor = {}
         
        ## calc coms for com-explosion and labeling
        coms = {}
        ligandsCOMS = {}
        chainsCOMS = {}
         
        for obj in selected:
            if typeOfExplosion == 'com':
                ## calculate com of obj
                coms[obj] = calc_COM(obj)
            s = s + ' '+ obj
        if not complex: 
            cmd.create('_all_obj', s)
            dim = cmd.get_extent('_all_obj')
            complexXYZ = calc_COM('_all_obj')
            if typeOfExplosion == 'com':
                ## calc com of complex from objects
                trans = calcTransFac('_all_obj')
            else:
                ## calculate translation vector
                transVec = transAxes('_all_obj') 
                print "Translation vector:", transVec
            cmd.frame(1)
            cmd.mview('store', object = 'all')
            cmd.delete('_all_obj')
                 
        else:
            dim = cmd.get_extent(complex)
            complexXYZ = calc_COM(complex)
            if typeOfExplosion == 'com':
                ## calc com of complex from objects
                trans = calcTransFac(complex)
                 
            else:
                ## calculate translation vector
                transVec = transAxes(complex) 
                print "Translation vector:",transVec
             
        for obj in selected:
            ##get chains of complex
            chains[obj] = cmd.get_chains(obj)
            for ch in chains[obj]:
                chainsAndObj[obj+ '_chain' + ch] = obj
                chainAndLigand[obj+ '_chain' + ch] = []
                 
            ## create objects for all chains and ligands
            if typeOfExplosion == 'com':
                cNames_new, chainAndLabel_new, chainAndLigand_new, \
                ligandAndChain_new, ligandsCOMS_new, chainsCOMS_new, f, \
                label_objects_new, objColor_new = \
                    create_objects(chains[obj], obj, storedLigands, chainsCOMS, 
                                    complexXYZ, dim, colorBinding, chainAndLigand)
                if ligandsCOMS_new:
                    ligandsCOMS.update(ligandsCOMS_new)
                if chainsCOMS_new:
                    chainsCOMS.update(chainsCOMS_new)
                     
            if typeOfExplosion == 'canonical':
                cNames_new, chainAndLabel_new, chainAndLigand_new, \
                ligandAndChain_new, f, label_objects_new, objColor_new=  \
                    create_objects(chains[obj], obj, storedLigands, chainsCOMS, 
                                    complexXYZ, dim, colorBinding, chainAndLigand, 
                                    typeOfExplosion='canonical')
                                     
            cNames = cNames + cNames_new
            if chainAndLabel_new:
                chainAndLabel.update(chainAndLabel_new)
            if chainAndLigand_new:
                chainAndLigand.update(chainAndLigand_new)
            if ligandAndChain_new:
                ligandAndChain.update(ligandAndChain_new)
            if label_objects_new:
                label_objects.update(label_objects_new)
            if objColor_new:
                objColor.update(objColor_new)
         
        ## color contact sites between chains
        if colorBinding == 'contact':
            contactColors = color_contact(cNames, objColor, colorBinding)   
        elif colorBinding == 'chain':
            color_contact(cNames, objColor, colorBinding)   
         
        ## color binding site
        if colorBinding != 'none':
            for cl in chainAndLigand.iteritems():
                for l in cl[1]:
                    if colorBinding == 'contact' and contactColors:
                        color_binding(cl[0], l[:-1], contactColors)
                    elif colorBinding == 'chain':
                        color_binding(cl[0], l[:-1])        
        ## set best view
        cmd.frame(f)
        view_objects = best_view_objects()
        best_view(view_objects, 'chain', '10')
        cmd.mview('store', object='all')    
         
        ## store chain objects for movie 
        store_view(group = True, all = True)
         
        ## labeling for canonical explosion
        if typeOfExplosion == 'canonical':
            chainslist = cNames
            ligandslist = ligandAndChain.keys()
            chainslist.extend(ligandslist)
            chainAndLabel, label_objects = \
                    calc_label_position_flush(chainslist, transVec, 1, 
                                                chainAndLabel)
            for ch in cNames:
                if ch in chainAndLigand.keys() and chainAndLigand[ch]:
                    for l in chainAndLigand[ch]:
                        cmd.group(l, label_objects[ch]+  " " + \
                        label_objects[l[:-1]] + " " + "label" + ch + " " + \
                        "label" + l[:-1] + " " + '_'+ch + '_label' + " " +\
                        '_'+l[:-1] + '_label', 'add') 
                cmd.group(ch + '_', "label"+ch + " "+'_'+ch + '_label','add')
            cmd.frame(1)
            cmd.zoom('all', complete=1)
            cmd.mview('store', object='all')
            if not showlabels:
                cmd.scene('on', 'store')
                cmd.mview('store', scene='on')
 
        # show labels
        if typeOfExplosion == 'com':
            cmd.frame(f-20)
            if not showlabels:
                show_labels(label_objects)
                cmd.scene('on', 'store')
                cmd.mview('store', scene='on')
                cmd.mview('store', object='all')
         
        f = f + 20
        cmd.frame(f)    
        cmd.zoom('all', complete=1)
        cmd.mview('store', object='all')
        if not showlabels:
            cmd.mview('store', scene='on')
         
        ## hide labels
        f = f + 1
        if not showlabels:
            cmd.frame(f)
            cmd.hide('labels')
            cmd.hide('dashes')
            cmd.scene('off', 'store')
            cmd.zoom('all',complete=1)
            cmd.mview('store', scene='off')
             
        f = f + 10     
        if len(selected) > 1:
            if typeOfExplosion == 'canonical':
                i = 1
            cmd.frame(f)
            ''' separate objects '''
            for obj in chains.keys():
                create_complex(chains, obj)
                     
                for obj2 in chains.keys():
                    if obj != obj2:
                        create_complex(chains, obj2)
                         
                        ## check if obj is colliding with any other obj
                        while isColliding('(_sel%s)'%obj, '(_sel%s)'%obj2):
                            if typeOfExplosion == 'com':
                                ## if collision, translate all chains in selection
                                for cname in cNames:
                                    com_translation(cname, chainAndLigand, 
                                                    complexXYZ, 
                                                    coms[chainsAndObj[cname]], 
                                                    trans, f)
                                 
                            else:
                                for o in chains.keys():
                                    for ch in chains[o]:
                                        ch = obj+ '_chain' + ch
                                        canonical_translation(ch, i, transVec, 
                                                    chainAndLigand, label_objects)
                                    i = i + 1
                            create_complex(chains, obj)
                            create_complex(chains, obj2)
                store_view(group = True, all = True)
                 
            f = f+ 30
            cmd.frame(f)
            store_view(group = True, all = True)
                 
            cmd.delete('(_sel%s)'%obj)
            cmd.delete('(_sel%s)'%obj2)
             
            print 'Preparation:', time.clock() - start_time, 'seconds'
         
        ''' translate objects individually '''
        for obj in chains.keys():
            if len(selected) > 1:
                f = f + 10
            sel1=''
            for c in chains[obj]:
                sel1 = sel1 + obj+'_chain' + c + ' '
            cmd.select('_'+obj, sel1)
            if typeOfExplosion == 'com':
                trans = calcTransFac('_'+obj)
            else:
                trans = transAxes('_'+obj)
             
            ## get object's chains and ligands
            objChains = [chain for chain, o in chainsAndObj.items() if o == obj]
            ligandAndChain_obj= {}
            for ch in objChains:
                if ch in ligandAndChain.values():
                    ligandAndChain_obj.update({ligand: c for ligand, 
                                        c in ligandAndChain.items() if c == ch})
             
            ## translation of object
            if typeOfExplosion == 'com':
                if not complex:
                    f = com_explosion('_'+obj, label_objects, objChains,
                                chainAndLigand, ligandAndChain_obj, chainsCOMS, 
                                ligandsCOMS, trans, com = coms[obj], frame = f)
                                 
                else:
                    f = com_explosion('_'+obj, label_objects, objChains,
                                chainAndLigand, ligandAndChain_obj, chainsCOMS, 
                                ligandsCOMS, trans, complex = complex,
                                com = complexXYZ, frame = f)
 
            else:
                if complex:
                    f = canonical_explosion('_'+obj, label_objects, objChains, 
                            chainAndLigand, trans,ligandAndChain_obj, complex, f)
                             
                else:
                    f = canonical_explosion('_'+obj, label_objects, objChains, 
                            chainAndLigand,trans, ligandAndChain_obj, frame = f)
                 
                 
            cmd.delete('_'+obj)
 
    if len(selected) > 1:
        f = f - 25
    cmd.frame(f)
    view_objects =best_view_objects(True)
    best_view(view_objects, 'chain', '10')
    cmd.zoom('all', complete=1)
    if not showlabels and len(selected) == 1:
    	#f = f -29
        cmd.frame(f)
        show_labels(label_objects)
        cmd.scene('on2', 'store')
        cmd.mview('store', scene='on2')
         
    print 'Explosion of', selected, time.clock() - start_time, 'seconds'
     
# video muss zweimal gespeichert werden (bei erstem mal noch wackeln in 
#           letztem gespeicherten frame) '''    
def relabel(selected, newLabels):
    '''DESCRIPTION:
        rename an selected object and its label by a new label
    '''
    obj = cmd.get_names("objects")
    x = ''
     
    selectedObjects = selected.split()
    labels = newLabels.split()
    ## relabel label
    for i in range(0, len(selectedObjects)):
        cmd.label("label" + selectedObjects[i], "'%s'" %labels[i])
        ## rename all objects
        for o in obj:
            if selectedObjects[i] in o:
                x = o
                x = x.replace(selectedObjects[i], labels[i])
                cmd.set_name(o, x)
                cmd.scene('on','update')
                cmd.scene('off','update')
                cmd.scene('on2','update')   
             
def reorient_explosion(frame=1):
    '''DESCRIPTION:
        if called, the orientation of the movie is set to the orientation of 
        given frame from frame till end. '''
    cmd.zoom('all', complete = 1)
    x = cmd.get_view(quiet=1)
    frames = cmd.count_frames()
    for f in range(int(frame),frames+1):
        cmd.frame(f)
        cmd.set_view(x)
        cmd.zoom('all', complete = 1)
        cmd.mview('store')
         
def renew_representation(selection, representation):
    '''DESCRIPTION:
        show selection in given representation for complete movie
    '''
    selection = selection.split()
    for f in range(0, cmd.count_frames()):
        cmd.frame(f)
        for s in selection:
            cmd.show_as(representation, s)
            labels_on()
            cmd.scene('on', 'update')
            cmd.scene('on2', 'update')
         
    hide_labels()
    cmd.scene('off', 'update')
             
def hide_labels():
    cmd.hide('labels')
    cmd.hide('dashes')
     
def labels_on():
    cmd.show('dashes')
    for obj in cmd.get_names('objects'):
        if obj.startswith('label'):
            cmd.show('labels', obj)
         
def remove_ligand(name):
    ''' DESCRIPTION:
        Remove specific solvent AFTER explosion
    '''
    for obj in cmd.get_names('objects'):
        if name in obj and obj[-1] != '_':
            cmd.delete(obj)
             
cmd.extend('explosion', explosion)
cmd.extend('relabel', relabel)
cmd.extend('reorient_explosion', reorient_explosion)
cmd.extend('renew_representation', renew_representation)
cmd.extend('show_labels', labels_on)
cmd.extend('hide_labels', hide_labels)
cmd.extend('remove_ligand', remove_ligand)