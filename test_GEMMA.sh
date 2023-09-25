./gemma --help


# Kinship
/common/jungj2/miniconda3/bin/gemma \
    -g example/mouse_hs1940.geno.txt.gz \
    -p example/mouse_hs1940.pheno.txt \
    -gk 1 \
    -o mouse_hs1940

# Perform Eigen-Decomposition of the Relatedness Matrix
/common/jungj2/miniconda3/bin/gemma \
    -g example/mouse_hs1940.geno.txt.gz \
    -p example/mouse_hs1940.pheno.txt \
    -k output/mouse_hs1940.cXX.txt \
    -eigen \
    -o mouse_hs1940

# LMM
/common/jungj2/miniconda3/bin/gemma \
    -g example/mouse_hs1940.geno.txt.gz \
    -p example/mouse_hs1940.pheno.txt \
    -a example/mouse_hs1940.anno.txt \
    -k output/mouse_hs1940.cXX.txt \
    -lmm 2 \
    -o mouse_hs1940_CD8_lmm


# mvLMM
/common/jungj2/miniconda3/bin/gemma \
    -g example/mouse_hs1940.geno.txt.gz \
    -p example/mouse_hs1940.pheno.txt \
    -a example/mouse_hs1940.anno.txt \
    -k output/mouse_hs1940.cXX.txt \
    -lmm \
    -n 1 6 \
    -o mouse_hs1940_CD8_mvlmm


#####################################################################
#####################################################################
#####################################################################

/common/jungj2/miniconda3/bin/gemma \
    -g /common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.geno.txt \
    -p /common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.pheno.txt \
    -gk 1 \
    -o mouse100

# Perform Eigen-Decomposition of the Relatedness Matrix
/common/jungj2/miniconda3/bin/gemma \
    -g /common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.geno.txt \
    -p /common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.pheno.txt \
    -k output/mouse100.cXX.txt \
    -eigen \
    -o mouse100

# mvLMM
/common/jungj2/miniconda3/bin/gemma \
    -g /common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.geno.txt \
    -p /common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.pheno.txt \
    -k /common/jungj2/miniconda3/envs/r_env/lib/R/library/gemma2/extdata/mouse100.cXX.txt \
    -lmm \
    -n 1 6 \
    -o mouse100


    