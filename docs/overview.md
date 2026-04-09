# Visão Geral do Sistema

## 1. Descrição

Este projeto implementa um simulador de sistemas gravitacionais N-body sob a mecânica clássica newtoniana, com foco em estabilidade numérica
eficiência computacional e visualização em tempo real.

O sistema permite a simulação de múltiplos corpos interagindo gravitacionalmente, com suporte a diferentes estratégias de cálculo de força e integração temporal.

---

## 2. Escopo e Premissas

O modelo físico adotado assume:
- Gravitação universal segundo a mecânica newtoniana
- Sistema isolado (sem forças externas)
- Corpos pontuais (sem volume ou colisão física explícita)
- Uso de softening gravitacional para evitar singularidades numéricas

Não são considerados:
- Efeitos relativísticos
- Interações não gravitacionais
- Modelos de colisão ou fusão entre corpos

---

## 3. Arquitetura Geral
O sistema é estruturado em três componentes principais:

### 3.1 Motor de Simulação (CPU)

Responsável por:

- Cálculo das forças gravitacionais
- Integração temporal do sistema
- Atualização do estado físico (posição, velocidade)

### 3.2 Sistema de Renderização (GPU)

Responsável por:

- Conversão de coordenadas de mundo para espaço de tela
- Renderização eficiente das trajetórias
- Execução de shaders procedurais para construção geométrica


## 6. Propriedades do Sistema

## 7. Limitações

As principais limitações incluem:

- Erro numérico acumulado devido à discretização temporal
- Dependência da escolha do timestep (dt)
- Aproximações introduzidas pelo algoritmo de Barnes–Hut
- Ausência de tratamento físico de colisões