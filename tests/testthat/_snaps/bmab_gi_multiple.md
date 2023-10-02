# bmab_gi_multiple() gives correct output

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["alpha", "beta", "gi", "Sigma", "n", "stage"]
        },
        "class": {
          "type": "character",
          "attributes": {},
          "value": ["data.frame"]
        },
        "row.names": {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 3]
        },
        "params": {
          "type": "list",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["alpha_start", "beta_start", "Sigma_start", "n_start", "gamma", "N", "num_actions", "tol"]
            }
          },
          "value": [
            {
              "type": "double",
              "attributes": {},
              "value": [2]
            },
            {
              "type": "double",
              "attributes": {},
              "value": [1]
            },
            {
              "type": "double",
              "attributes": {},
              "value": [2]
            },
            {
              "type": "double",
              "attributes": {},
              "value": [3]
            },
            {
              "type": "double",
              "attributes": {},
              "value": [0.5]
            },
            {
              "type": "double",
              "attributes": {},
              "value": [20]
            },
            {
              "type": "double",
              "attributes": {},
              "value": [2]
            },
            {
              "type": "double",
              "attributes": {},
              "value": [0.0005]
            }
          ]
        },
        "gi_matrix": {
          "type": "double",
          "attributes": {
            "dim": {
              "type": "integer",
              "attributes": {},
              "value": [2, 2]
            },
            "dimnames": {
              "type": "list",
              "attributes": {
                "names": {
                  "type": "character",
                  "attributes": {},
                  "value": ["beta", "alpha"]
                }
              },
              "value": [
                {
                  "type": "character",
                  "attributes": {},
                  "value": ["1", "2"]
                },
                {
                  "type": "character",
                  "attributes": {},
                  "value": ["2", "3"]
                }
              ]
            }
          },
          "value": [0.70600586, 0.53586272, 0.77735112, "NA"]
        },
        "gi_matrix_ns": {
          "type": "double",
          "attributes": {
            "dim": {
              "type": "integer",
              "attributes": {},
              "value": [2, 2]
            },
            "dimnames": {
              "type": "list",
              "attributes": {
                "names": {
                  "type": "character",
                  "attributes": {},
                  "value": ["n", "Sigma"]
                }
              },
              "value": [
                {
                  "type": "character",
                  "attributes": {},
                  "value": ["3", "4"]
                },
                {
                  "type": "character",
                  "attributes": {},
                  "value": ["2", "3"]
                }
              ]
            }
          },
          "value": [0.70600586, 0.53586272, "NA", 0.77735112]
        }
      },
      "value": [
        {
          "type": "integer",
          "attributes": {},
          "value": [2, 2, 3]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [1, 2, 1]
        },
        {
          "type": "double",
          "attributes": {},
          "value": [0.70600586, 0.53586272, 0.77735112]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [2, 2, 3]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [3, 4, 4]
        },
        {
          "type": "integer",
          "attributes": {},
          "value": [0, 1, 1]
        }
      ]
    }

