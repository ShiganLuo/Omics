Traceback (most recent call last):
  File "/data/pub/zhousha/env/mutation_0.1/1d458d8737941359359538b86c40e1d5_/bin/umi_tools", line 11, in <module>
    sys.exit(main())
  File "/data/pub/zhousha/env/mutation_0.1/1d458d8737941359359538b86c40e1d5_/lib/python3.6/site-packages/umi_tools/umi_tools.py", line 61, in main
    module.main(sys.argv)
  File "/data/pub/zhousha/env/mutation_0.1/1d458d8737941359359538b86c40e1d5_/lib/python3.6/site-packages/umi_tools/dedup.py", line 293, in main
    threshold=options.threshold)
  File "/data/pub/zhousha/env/mutation_0.1/1d458d8737941359359538b86c40e1d5_/lib/python3.6/site-packages/umi_tools/network.py", line 403, in __call__
    clusters = self.UMIClusterer(counts, threshold)
  File "/data/pub/zhousha/env/mutation_0.1/1d458d8737941359359538b86c40e1d5_/lib/python3.6/site-packages/umi_tools/network.py", line 368, in __call__
    min(len_umis), max(len_umis)))
AssertionError: not all umis are the same length(!):  19 - 20