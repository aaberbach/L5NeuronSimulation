for i in 0.45 0.55 0.65 0.75 0.85 0.95 2.4 2.5 2.6 2.7 2.8 2.9 3 #0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3
do
    python build_network.py $i
    python run_network.py $i
done