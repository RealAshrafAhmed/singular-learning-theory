def __set_param(probs, weights, name, index, value):
    m_probs = probs.copy()
    m_weights = weights.copy()
    if name == "probs":
        probs[index] = value
    elif name == "weights":
        weights[index] = value
    else:
        raise Exception(f"Not support name={name}")

    # print(f"changed probs={probs}, weights={weights} to probs={m_probs}, weights={m_weights} using {name}, {index}, {value}")
    return probs, weights

def __set_params(probs, weights, params):
    m_probs = probs.copy()
    m_weights = weights.copy()
    for param in params:
        name = param["name"]
        index = param["index"]
        value = param["value"]
        m_probs, m_weights = __set_param(m_probs, weights, name, index, value)

    return m_probs, m_weights