context("Write and read probs")

test_that("write_probs and read_probs work", {

    if(isnt_karl()) skip("this test only run locally")

    library(qtl2geno)
    iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
    iron <- iron[,c(18,19,"X")]
    map <- insert_pseudomarkers(iron, step=0.5)
    pr <- calc_genoprob(iron, map, error_prob=0.002)

    db <- "test_iron_probs"

    write_probs(db, pr, map)

    pr18 <- read_probs(db, chr="18")
    expect_equal(pr18, pr[,18])

    prX <- read_probs(db, chr="X", pos=c(29.8, 32.8))
    expected <- pr[,"X"]
    expected$X <- expected$X[,,2:7]
    expect_equal(prX, expected)

})

test_that("write_probs and read_probs work for backcross", {

    if(isnt_karl()) skip("this test only run locally")

    library(qtl2geno)
    grav <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
    grav <- grav[,4:5]
    map <- insert_pseudomarkers(grav, step=1, stepwidth="max")
    pr <- calc_genoprob(grav, map, error_prob=0.002)

    db <- "test_grav_probs"

    write_probs(db, pr, map)

    pr4 <- read_probs(db, chr="4")
    expect_equal(pr4, pr[,4])

    pr5 <- read_probs(db, chr="5", pos=c(29.8, 32.8))
    expected <- pr[,"5"]
    expected[["5"]] <- expected[["5"]][,,40:42]
    expect_equal(pr5, expected)

})
