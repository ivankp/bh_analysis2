#!/bin/bash

for P in I
do

rootmv H2j_${P}_sherpa.root:MUR0.5_MUF0.5_PDF11000_JetAntiKt4 \
       H2j_${P}_sherpa.root:FacHt4_RenHt4_PDFCT10nlo_JetAntiKt4

rootmv H2j_${P}_sherpa.root:MUR0.5_MUF1_PDF11000_JetAntiKt4 \
       H2j_${P}_sherpa.root:FacHt2_RenHt4_PDFCT10nlo_JetAntiKt4

rootmv H2j_${P}_sherpa.root:MUR1_MUF0.5_PDF11000_JetAntiKt4 \
       H2j_${P}_sherpa.root:FacHt4_RenHt2_PDFCT10nlo_JetAntiKt4

rootmv H2j_${P}_sherpa.root:MUR1_MUF1_PDF11000_JetAntiKt4 \
       H2j_${P}_sherpa.root:FacHt2_RenHt2_PDFCT10nlo_JetAntiKt4

rootmv H2j_${P}_sherpa.root:MUR1_MUF2_PDF11000_JetAntiKt4 \
       H2j_${P}_sherpa.root:FacHt_RenHt2_PDFCT10nlo_JetAntiKt4

rootmv H2j_${P}_sherpa.root:MUR2_MUF1_PDF11000_JetAntiKt4 \
       H2j_${P}_sherpa.root:FacHt2_RenHt_PDFCT10nlo_JetAntiKt4

rootmv H2j_${P}_sherpa.root:MUR2_MUF2_PDF11000_JetAntiKt4 \
       H2j_${P}_sherpa.root:FacHt_RenHt_PDFCT10nlo_JetAntiKt4

done
