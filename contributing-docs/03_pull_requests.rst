Pull Request Guidelines
========================

Ce document explique comment créer des Pull Requests et détaille les standards de code attendus lors de leur mise en œuvre.

Avant de soumettre
------------------

Checklist pré-soumission
~~~~~~~~~~~~~~~~~~~~~~~~

Avant de soumettre une PR depuis ton fork, vérifie qu'elle respecte les points suivants :

✔️ **Tests requis**

Inclure des tests (doctests, tests unitaires avec pytest, ou les deux).
Les tests doivent couvrir les nouvelles fonctionnalités et les cas limites.

✔️ **Build ReadTheDocs vert**

Le build de la documentation doit passer sans erreur.
Les mainteneurs ne fusionneront jamais une PR qui casse le lint ou la documentation.

✔️ **Conversations résolues**

Toutes les discussions doivent être résolues avant que la PR soit fusionnée.

✔️ **Rebase plutôt que merge**

Il est demandé de rebaser fréquemment ta PR pour garder un historique propre et faciliter la revue.
Tous les conflits doivent être résolus.

✔️ **Fusion en "Squash and Merge"**

Peu importe le nombre de commits pendant la review, la PR sera fusionnée en un seul commit.
Les mainteneurs peuvent te demander de nettoyer ou regrouper les commits avant la fusion.

✔️ **Licence MIT obligatoire**

Tout nouveau fichier doit commencer par l'en-tête de licence MIT.

✔️ **Code + tests + docs dans la même PR**

Si tu ajoutes une fonctionnalité, les documents doivent être mis à jour dans la même PR :

- Docstrings dans le code
- Documentation Sphinx si nécessaire
- Fichiers README si applicable

✔️ **PR petites et focalisées**

Ne mélange pas refactorings et nouvelles features.
Les petites PR se review beaucoup mieux et facilitent le cherry-pick pour les releases correctives.

Pour les grosses modifications :

1. Créer un Draft global pour discussion
2. Envoyer ensuite plusieurs petites PR dérivées

✔️ **Exécuter les tests localement**

Les tests suivent la même arborescence que le code.

Exemple : changements dans ``spatialetl-core/spatialetl/coverage/`` ⇒ tests dans ``spatialetl-core/spatialetl/coverage/tests/``

Lancer les tests avec UV :

.. code-block:: bash

   # Tests d'un module spécifique
   uv run pytest spatialetl-core/spatialetl/coverage/tests/

   # Tous les tests
   uv run pytest

✔️ **Tester sur Python 3.9**

Version minimale supportée : Python 3.9.
Certaines fonctionnalités récentes (match/case, nouveaux types) ne fonctionnent pas sur cette version.

✔️ **Messages de commit conformes**

Format recommandé : ``[Type] Description courte``

Types acceptés :

- ``[Feat]`` : Nouvelle fonctionnalité
- ``[Fix]`` : Correction de bug
- ``[Docs]`` : Documentation uniquement
- ``[Chore]`` : Maintenance, dépendances
- ``[Refactor]`` : Refactoring sans changement fonctionnel
- ``[Test]`` : Ajout ou modification de tests

Exemple : ``[Feat] add AROME forecast data provider``

Processus de review
-------------------

Résolution des conversations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Une PR est mergeable uniquement lorsque **toutes les conversations sont résolues**.

Cela permet :

- Une visibilité claire de l'état de la PR
- Une boucle review/merge plus rapide
- De limiter l'usage de "Request changes" aux cas vraiment bloquants

Ce que les reviewers attendent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Code lisible et bien documenté
- Tests pertinents et qui passent
- Documentation à jour

Après la review
~~~~~~~~~~~~~~~

Si des modifications sont demandées :

1. Apporter les corrections
2. Commit et push sur ta branche
3. Répondre aux commentaires pour indiquer que c'est fait
4. Marquer les conversations comme résolues si approprié

Les mainteneurs fusionneront dès que tout est validé.