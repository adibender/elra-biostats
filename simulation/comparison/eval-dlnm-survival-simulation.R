library(batchtools)
reg <- loadRegistry("dlnm-tv-surv-registry/")
library(ggplot2)
theme_set(theme_bw())

### eval dlnm_ped
id_pam_ped <- findExperiments(prob.name="dlnm_sim_ped", algo.name="pam_dlnm_ped")
res_pam_ped <- reduceResultsDataTable(ids=findDone(id_pam_ped[,1])) %>%
	as_tibble() %>%
	tidyr::unnest()

av_pam_ped <- res_pam_ped %>%
	group_by(x, lag2) %>%
	summarize_at(vars(fit, truth), mean ) %>%
	arrange(x, lag2) %>%
	tidyr::gather(type, value, fit:truth)

jpeg("dlnm_ped.jpg", width=700, height=300)
ggplot(av_pam_ped, aes(x=x, y=lag2, fill=value)) +
	geom_tile() +
	geom_contour(aes(z=value), col="grey70") +
	scale_fill_gradient2(low="firebrick", high="steelblue") +
	facet_wrap(~type)
dev.off()


#### PED PAM DLNM TV
id_pam_dlnm_ped_tv <- findExperiments(
	prob.name="sim_dlnm_ped_tv",
	algo.name="pam_dlnm_ped_tv")
res_pam_dlnm_ped_tv <- reduceResultsDataTable(ids=findDone(id_pam_dlnm_ped_tv[,1])) %>%
	as_tibble() %>%
	tidyr::unnest()

av_tv <- res_pam_dlnm_ped_tv %>%
	group_by(x, lag2, time_df) %>%
	summarize_at(vars(fit, truth), mean ) %>%
	arrange(time_df)

ggsave("simulation/av_tv_dlnm_pam.pdf", width=6, height=10)

library(tweenr)
library(gganimate)
anim_df <- av_tv %>%
	mutate(
		group = row_number(),
		time = round(time_df),
		ease = "cubic-in-out") %>%
	select(-time_df) %>%
	tidyr::gather(type, value, fit:truth) %>%
	arrange(time)
# tween_3d <- tween_elements(anim_df, time="time_df", group="group", ease="ease",
# 	nframes=400)
gg_3d <- ggplot(anim_df, aes(x=x, y=lag2, fill=value, z=value, frame=time)) +
	geom_tile() +
	geom_contour(aes(group=time), col="grey70") +
	scale_fill_gradient2(low="firebrick", high="steelblue") +
	facet_wrap(~type)
gganimate(gg_3d, "fit-vs-truth-ped.gif", interval=0.3, ani.width=700, height=300)
