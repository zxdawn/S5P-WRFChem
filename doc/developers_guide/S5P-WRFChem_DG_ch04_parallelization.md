<!-- BEGIN COMMENT -->

 [<< Previous Chapter](S5P-WRFChem_DG_ch03_coding.md)- [Home](README.md)

<!-- END COMMENT -->

# 4. Parallelization

The S5P-WRFChem uses [Dask](https://dask.org/) to parallelize the retrieval.

Two parameters in the `setting.txt` controls it:

```
n_workers = 4                      ; Dask: number of workers to start
threads_per_worker = 1             ; Dask: number of threads per each worker
```

As Xin Zhang is the dask beginner, if you have any good idea, PRs are welcome ;)

<!-- BEGIN COMMENT -->

 [<< Previous Chapter](S5P-WRFChem_DG_ch03_code)- [Home](README.md)<br>
S5P-WRFChem User's Guide (c) 2021<br>

<!-- END COMMENT -->