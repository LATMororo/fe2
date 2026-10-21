from _fe2_material import MaterialFactory
from pyprnn.prnn   import PRNN
from pyprnn.utils  import StressStrainDataset, Trainer

import numpy as np
import torch

J2 = '''
  material = {
    type        = "J2"; 
    rank        = 3; 
    state       = "Solid";
    young       = 79500.; 
    poisson     = 0.33; 
    yield       = "90.0740 * ( 1. - exp(-113.0268 * x)) + 444.1184";
    rmTolerance = 1.e-12; 
    rmMaxIter   = 100;
  };'''

torch.set_default_dtype(torch.float64)
#torch.manual_seed(42)

# Load dataset

dataset = StressStrainDataset(
        'monotonic_loading.out', list(range(6)), list(range(6,12)), seq_length=61)

# Split data into training and validation ones

tset, vset = torch.utils.data.random_split(
        dataset, [0.3, 0.7], generator=torch.Generator().manual_seed(42))

tloader = torch.utils.data.DataLoader(tset, batch_size=3, shuffle=True)
vloader = torch.utils.data.DataLoader(vset, batch_size=10, shuffle=False)

# Initialize prnn model

factory = MaterialFactory(J2)

prnn = PRNN(n_features=6, n_outputs=6, n_matpts=2, material_factory=factory).double() # double precision 

# Train network
trainer = Trainer(prnn, optimizer=torch.optim.Adam(prnn.parameters(), lr=1e-3))

trainer.train(tloader, tloader, epochs=5, patience=650)

trainer.save('prnn_J2_material.pth')
