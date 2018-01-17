library(transport)

n <- 5000
site <- rbinom(n, 1, .5)

race <- rbinom(n,1, .4 + (.2 * site))
crime <- rnorm(n, .1 * site , 1)
discrimination <- rnorm(n, 1 + (.2 * site ), 1)

# Instrument
voucher <- rbinom(n, 1, .5)

# Exposure
move0 <- rbinom(n, 1, plogis(-log(1.6) -
                               log(1.1) * crime - log(1.3) * discrimination))
move1 <- rbinom(n, 1, plogis(-log(1.6) + log(4) - log(1.1) * crime -
                               log(1.3) * discrimination))
move <- ifelse(voucher == 1, move1, move0)

# Outcomes
inschoola0 <- rbinom(n, 1, plogis(log(1.6) +
                           (log(1.9) * move0) - log(1.3) * discrimination -
                            log(1.2) * race + log(1.2) * race * move0))
inschoola1 <- rbinom(n, 1, plogis(log(1.6) +
                           (log(1.9) * move1) - log(1.3) * discrimination -
                            log(1.2) * race + log(1.2) * race * move1))
inschoola <- ifelse(voucher == 1, inschoola1, inschoola0)

dat <- data.frame(w2 = crime , w3 = discrimination, w1 = race , site = site ,
                  a = voucher, z = move, y = inschoola)
wmat <- dat[, c("w1", "w2", "w3")]


amodel <- "a ~ 1"
sitemodel <- "site ~ w1 + w2 + w3 "
zmodel <- "z ~ a + w2 + w3"
outmodel <- "y ~ z + w1 + w3 + z:w1"
outmodelnoz <- "y ~ a + w1 + w3 + a:w1"
q2model <- "w1 + w2 + w3 "

# TODO: provide s_awz_model formula?
ittate_est <- ittatetmle(a = dat$a, z = dat$z, y = dat$y, site = dat$site,
                         w = wmat, aamodel = amodel, asitemodel = sitemodel,
                         #azmodel = zmodel,
                         s_awz_model = zmodel,
                         aoutmodel = outmodel,
                         aq2model = q2model)

cate_est <- catetmle(ca = dat$a, cz = dat$z, cy = dat$y, csite = dat$site,
                     cw = wmat, csitemodel = sitemodel, czmodel = zmodel,
                     coutmodel = outmodel, cq2model = q2model)

eate_est <- eatetmle(a = dat$a, z = dat$z, y = dat$y, site = dat$site, w = wmat,
                     nsitemodel = sitemodel, nzmodel = zmodel,
                     noutmodel = outmodel)